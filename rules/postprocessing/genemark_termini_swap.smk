"""
Swap BRAKER CDS termini with GeneMark's when intron chains agree.

When BRAKER and GeneMark independently predict the same CDS intron chain
for a transcript but disagree on start/stop positions, GeneMark's termini
are correct ~90% of the time (benchmarked on A. thaliana ET/EP/ETP:
3,245-4,008 transcripts improved, only 139-265 worsened per mode).

The intron chain agreement between two independent gene finders is the
confidence signal.  Only multi-exon transcripts are swapped (single-exon
matches are vacuous).
"""


def _termini_swap_genemark_input(wildcards):
    """Return the GeneMark GTF path for the sample's pipeline mode."""
    sample = wildcards.sample
    mode = get_braker_mode(sample)
    if mode == 'es':
        return f"output/{sample}/GeneMark-ES/genemark.gtf"
    elif mode == 'ep':
        return f"output/{sample}/GeneMark-EP/genemark.gtf"
    elif mode == 'et':
        return f"output/{sample}/genemark/genemark.gtf"
    else:  # etp, isoseq, dual
        return f"output/{sample}/GeneMark-ETP/genemark.gtf"


rule genemark_termini_swap:
    input:
        braker_gtf="output/{sample}/braker.tsebra.gtf",
        genemark_gtf=_termini_swap_genemark_input,
    output:
        swapped_gtf="output/{sample}/braker.termini_swapped.gtf",
        swap_log="output/{sample}/genemark_termini_swap.log",
    benchmark:
        "benchmarks/{sample}/genemark_termini_swap/genemark_termini_swap.txt"
    params:
        script=os.path.join(script_dir, "genemark_termini_swap.py"),
    threads: 1
    resources:
        mem_mb=int(config['slurm_args']['mem_of_node']) // int(config['slurm_args']['cpus_per_task']),
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        set -euo pipefail
        python3 {params.script} \
            --braker {input.braker_gtf} \
            --genemark {input.genemark_gtf} \
            --output {output.swapped_gtf} \
            2>&1 | tee {output.swap_log}
        """
