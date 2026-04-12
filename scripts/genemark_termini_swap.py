#!/usr/bin/env python3
"""Swap BRAKER CDS termini with GeneMark's when intron chains agree.

When BRAKER and GeneMark independently predict the same CDS intron chain
for a multi-exon transcript but disagree on start/stop codon positions,
GeneMark's termini are correct ~90% of the time (benchmarked on
A. thaliana ET/EP/ETP: 3,245-4,008 transcripts improved, only 139-265
worsened per mode).

Only multi-exon transcripts are considered -- for single-exon genes the
intron chain is empty and the match is vacuous.
"""

import argparse
import re
import sys
from collections import defaultdict


def parse_gtf_transcripts(gtf_path):
    """Parse a GTF into per-transcript CDS and stop_codon features.

    Returns:
        dict: transcript_id -> {
            'chrom': str,
            'strand': str,
            'cds': [(start, end), ...],        # sorted by start
            'stop_codons': [(start, end), ...],
            'lines': [(line_index, feature_type, raw_line), ...],
        }
    """
    transcripts = defaultdict(lambda: {
        'chrom': None, 'strand': None,
        'cds': [], 'stop_codons': [], 'lines': [],
    })

    with open(gtf_path) as fh:
        for i, line in enumerate(fh):
            if line.startswith('#') or line.strip() == '':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature not in ('CDS', 'stop_codon'):
                continue

            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]

            # Extract transcript_id
            m = re.search(r'transcript_id\s+"([^"]+)"', cols[8])
            if not m:
                continue
            tid = m.group(1)

            tx = transcripts[tid]
            tx['chrom'] = chrom
            tx['strand'] = strand
            if feature == 'CDS':
                tx['cds'].append((start, end))
            elif feature == 'stop_codon':
                tx['stop_codons'].append((start, end))

    # Sort CDS exons by start position
    for tid, tx in transcripts.items():
        tx['cds'].sort()
        tx['stop_codons'].sort()

    return transcripts


def merge_stop_into_cds(cds_list, stop_codons, strand):
    """Merge stop_codon features into adjacent CDS (AUGUSTUS convention).

    AUGUSTUS puts stop codons as separate features outside the last CDS.
    GeneMark includes the stop codon in the last CDS. This normalizes
    both to the GeneMark convention (stop included in CDS) so intron
    chain comparison and boundary comparison are consistent.

    Returns a new sorted CDS list with the stop codon merged.
    """
    if not stop_codons or not cds_list:
        return list(cds_list)

    merged = [list(c) for c in cds_list]
    merged.sort()

    for sc_start, sc_end in stop_codons:
        if strand == '+':
            # Stop codon should be at or just after the last CDS exon
            last = merged[-1]
            if sc_start == last[1] + 1:
                last[1] = sc_end
            elif sc_start <= last[1] and sc_end >= last[1]:
                last[1] = max(last[1], sc_end)
        else:
            # Minus strand: stop codon is at or just before the first CDS exon
            first = merged[0]
            if sc_end == first[0] - 1:
                first[0] = sc_start
            elif sc_end >= first[0] and sc_start <= first[0]:
                first[0] = min(first[0], sc_start)

    return [tuple(c) for c in merged]


def intron_chain(cds_list):
    """Compute the intron chain from sorted CDS exons.

    Returns tuple of (intron_start, intron_end) pairs.
    """
    if len(cds_list) < 2:
        return ()
    introns = []
    for i in range(len(cds_list) - 1):
        intron_start = cds_list[i][1] + 1
        intron_end = cds_list[i + 1][0] - 1
        introns.append((intron_start, intron_end))
    return tuple(introns)


def build_genemark_lookup(gm_transcripts):
    """Build lookup: (chrom, strand, intron_chain) -> (first_cds_start, last_cds_end).

    Only multi-exon transcripts are included.
    """
    lookup = {}
    for tid, tx in gm_transcripts.items():
        cds = merge_stop_into_cds(tx['cds'], tx['stop_codons'], tx['strand'])
        if len(cds) < 2:
            continue
        ic = intron_chain(cds)
        if not ic:
            continue
        key = (tx['chrom'], tx['strand'], ic)
        # First one wins (they almost always agree on termini)
        if key not in lookup:
            lookup[key] = (cds[0][0], cds[-1][1])
    return lookup


def swap_termini(braker_gtf, genemark_gtf, output_gtf):
    """Read braker GTF, swap termini where intron chains match GeneMark, write output."""

    # Parse GeneMark transcripts and build lookup
    gm_transcripts = parse_gtf_transcripts(genemark_gtf)
    gm_lookup = build_genemark_lookup(gm_transcripts)

    # Parse BRAKER transcripts
    br_transcripts = parse_gtf_transcripts(braker_gtf)

    # Determine which transcripts to swap and their new boundaries
    swap_info = {}  # tid -> (new_first_start, new_last_end)

    n_matched = 0
    n_swapped = 0
    n_unchanged = 0
    n_frame_skip = 0

    for tid, tx in br_transcripts.items():
        cds = merge_stop_into_cds(tx['cds'], tx['stop_codons'], tx['strand'])
        if len(cds) < 2:
            continue
        ic = intron_chain(cds)
        if not ic:
            continue
        key = (tx['chrom'], tx['strand'], ic)
        if key in gm_lookup:
            n_matched += 1
            gm_start, gm_end = gm_lookup[key]
            br_start = cds[0][0]
            br_end = cds[-1][1]
            if br_start != gm_start or br_end != gm_end:
                # Guard: if the start shift is not a multiple of 3, swapping
                # would change the reading frame through every downstream CDS
                # exon.  This means the two gene finders predicted different
                # ORFs that happen to share splice sites -- too risky to swap.
                if (gm_start - br_start) % 3 != 0:
                    n_frame_skip += 1
                    continue
                swap_info[tid] = (gm_start, gm_end)
                n_swapped += 1
            else:
                n_unchanged += 1

    print(f"[INFO] GeneMark lookup: {len(gm_lookup)} unique multi-exon intron chains",
          file=sys.stderr)
    print(f"[INFO] BRAKER multi-exon transcripts with intron chain match: {n_matched}",
          file=sys.stderr)
    print(f"[INFO] Transcripts swapped (termini differed): {n_swapped}", file=sys.stderr)
    print(f"[INFO] Transcripts unchanged (termini already agreed): {n_unchanged}",
          file=sys.stderr)
    if n_frame_skip:
        print(f"[WARN] Transcripts skipped (start shift not multiple of 3): {n_frame_skip}",
              file=sys.stderr)

    # Second pass: rewrite the GTF with swapped boundaries
    with open(braker_gtf) as fin, open(output_gtf, 'w') as fout:
        for line in fin:
            if line.startswith('#') or line.strip() == '':
                fout.write(line)
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                fout.write(line)
                continue

            feature = cols[2]
            if feature not in ('CDS', 'start_codon', 'stop_codon'):
                fout.write(line)
                continue

            m = re.search(r'transcript_id\s+"([^"]+)"', cols[8])
            if not m or m.group(1) not in swap_info:
                fout.write(line)
                continue

            tid = m.group(1)
            new_start, new_end = swap_info[tid]
            tx = br_transcripts[tid]
            strand = tx['strand']
            orig_cds = sorted(tx['cds'])

            if len(orig_cds) < 2:
                fout.write(line)
                continue

            first_cds = orig_cds[0]
            last_cds = orig_cds[-1]
            feat_start = int(cols[3])
            feat_end = int(cols[4])

            if feature == 'CDS':
                # Modify first CDS exon start
                if feat_start == first_cds[0] and feat_end == first_cds[1]:
                    cols[3] = str(new_start)
                    # Frame of the first CDS exon resets to 0 when start moves
                    if new_start != first_cds[0]:
                        cols[7] = '0'
                # Modify last CDS exon end
                elif feat_start == last_cds[0] and feat_end == last_cds[1]:
                    cols[4] = str(new_end)
                # Edge case: single-exon transcript that somehow matched
                # (shouldn't happen due to len>=2 guard, but be safe)
                elif (feat_start == first_cds[0] and feat_end == last_cds[1]
                      and first_cds == last_cds):
                    cols[3] = str(new_start)
                    cols[4] = str(new_end)
                    if new_start != first_cds[0]:
                        cols[7] = '0'

            elif feature == 'start_codon':
                if strand == '+':
                    cols[3] = str(new_start)
                    cols[4] = str(new_start + 2)
                else:
                    cols[3] = str(new_end - 2)
                    cols[4] = str(new_end)

            elif feature == 'stop_codon':
                if strand == '+':
                    # Stop codon is 3 bp after the last CDS end
                    # (in AUGUSTUS convention stop_codon is separate)
                    cols[3] = str(new_end + 1)
                    cols[4] = str(new_end + 3)
                else:
                    cols[3] = str(new_start - 3)
                    cols[4] = str(new_start - 1)

            fout.write('\t'.join(cols) + '\n')

    return n_swapped, n_matched, n_unchanged


def main():
    parser = argparse.ArgumentParser(
        description='Swap BRAKER CDS termini with GeneMark when intron chains agree.')
    parser.add_argument('--braker', required=True,
                        help='BRAKER GTF (from best_by_compleasm)')
    parser.add_argument('--genemark', required=True,
                        help='GeneMark GTF')
    parser.add_argument('--output', required=True,
                        help='Output GTF with swapped termini')
    args = parser.parse_args()

    n_swapped, n_matched, n_unchanged = swap_termini(
        args.braker, args.genemark, args.output)

    print(f"[DONE] Wrote {args.output}: {n_swapped} transcripts swapped, "
          f"{n_unchanged} already agreed, "
          f"{n_matched} total intron chain matches.", file=sys.stderr)


if __name__ == '__main__':
    main()
