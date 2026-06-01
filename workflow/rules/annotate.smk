import os
import sys
from shutil import rmtree
import glob
from math import ceil


include: "common.smk"


configfile: os.path.join(workflow.basedir, "../../config/config.yaml")

wildcard_constraints:
    sample="[^/]+",  # sample cannot contain forward slashes
    batch="stdin\.part_[0-9]+",

localrules:
    all,


if not os.path.exists("logs"):
    os.makedirs("logs")


def count_contigs_with_minlen(fasta_path, minlen):
    n = 0
    length = 0
    in_seq = False
    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if in_seq and length >= minlen:
                    n += 1
                in_seq = True
                length = 0
            else:
                length += len(line.strip())
        if in_seq and length >= minlen:
            n += 1
    return n


def calculate_nbatches(assembly_path, minlen, nseqs, max_batches=20):
    """Calculate number of batches needed for splitting assembly"""
    if not os.path.exists(assembly_path):
        raise ValueError(f"Assembly not found: {assembly_path}")
    n_filtered = count_contigs_with_minlen(assembly_path, minlen)
    if n_filtered <= 0:
        nbatches = 1
    else:
        nbatches = min(ceil(n_filtered / nseqs), n_filtered, max_batches)
    print(f"Detected {n_filtered} contigs >= {minlen} bp; splitting into {nbatches} parts")
    return nbatches


def get_batches(wildcards):
    """Return batch name list for a sample after checkpoint completes."""
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    return sorted([os.path.basename(b).replace(".fasta", "") for b in batches])


SAMPLES = config["sample"] if isinstance(config["sample"], list) else [config["sample"]]

abricates = expand(
    "{sample}/abricate/{tool}.tab",
    sample=SAMPLES,
    tool=[
        "argannot",
        "card",
        "ecoh",
        "ecoli_vf",
        "megares",
        "ncbi",
        "plasmidfinder",
        "resfinder",
        "vfdb",
    ],
)

outputs = (
    expand("{sample}/antismash/{sample}_antismash.gbk", sample=SAMPLES)
    + expand("{sample}/antismash/{sample}_antismash.tab", sample=SAMPLES)
    + expand("{sample}/amrfinder/{sample}_amrfinder.tab", sample=SAMPLES)
    + expand("{sample}/cazi/{sample}_cazi_overview.txt", sample=SAMPLES)
    + expand("{sample}/cazi/{sample}_cazi_substrate.out", sample=SAMPLES)
    + expand("{sample}/cazi/{sample}_annotated_cazymes_RPM.tsv", sample=SAMPLES)
    + abricates
)

if config.get("check_contigs", False):
    outputs += expand("{sample}/annotation/{sample}_metaerg.gff", sample=SAMPLES)


rule all:
    input:
        outputs,


checkpoint split_assembly:
    """Split assembly into chunks for parallel processing.

    Uses seqkit shuffle (--two-pass requires a seekable file, i.e. a real
    path on disk – not stdin) then filters by minimum length and splits into
    N parts.  The number of parts is pre-computed in Python so that the
    wildcard constraint on {batch} stays valid.
    """
    input:
        assembly=lambda wc: (
            config["assembly"][wc.sample]
            if isinstance(config["assembly"], dict)
            else config["assembly"]
        ),
    output:
        splitdir=directory("tmp/{sample}"),
        done=temp("tmp/{sample}/.split_done"),
    params:
        outdir=lambda wc: os.path.abspath("tmp/" + wc.sample),
        minlen=config["contig_annotation_thresh"],
        nbatches=lambda wc, input: calculate_nbatches(
            str(input.assembly),
            config["contig_annotation_thresh"],
            config.get("chunk_size_contigs", 200),
            config.get("max_batches", 20),
        ),
    container:
        "docker://pegi3s/seqkit:2.3.0"
    threads: 4
    resources:
        mem_mb=8000,
    log:
        e="{sample}/logs/split_assembly.e",
        o="{sample}/logs/split_assembly.o",
    shell:
        # NOTE: --two-pass requires a real file path (not stdin); input.assembly
        # must be a seekable file, which it always is when supplied via config.
        """
        set -euo pipefail
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        seqkit shuffle {input.assembly} --two-pass \
            | seqkit seq --min-len {params.minlen} --threads {threads} \
            | seqkit split2 \
                --by-part {params.nbatches} \
                --out-dir {params.outdir} \
                --force \
                > "{log.o}" 2>> "{log.e}"
        touch {output.done}
        """


rule annotate_orfs:
    """Run MetaERG on one assembly chunk.

    MetaERG's output_report.pl is known to fail in some versions
    (see https://github.com/xiaoli-dong/metaerg/issues/12).  We therefore
    accept either master.gff.txt (full success) or all.gff (partial success)
    and fail loudly only when neither exists.
    """
    container:
        config["docker_metaerg"]
    input:
        assembly="tmp/{sample}/{batch}.fasta",
    output:
        gff=temp("{sample}/annotation/{batch}/data/either_all_or_master.gff"),
        ffn=temp("{sample}/annotation/{batch}/data/cds.ffn"),
        faa=temp("{sample}/annotation/{batch}/data/cds.faa"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8 * 1024,
        runtime=lambda wildcards, attempt: attempt * 45,
    threads: 4
    params:
        metaerg_db_dir=config["metaerg_db_dir"],
        outdir=lambda wc: f"{wc.sample}/annotation/{wc.batch}",
    shell:
        # FIX 5: Simplified, non-redundant GFF fallback logic.
        # master.gff.txt  -> full success
        # data/all.gff    -> partial success (output_report.pl failed)
        # neither         -> hard failure
        """
        set +e
        metaerg.pl \
            --cpus {threads} \
            --dbdir {params.metaerg_db_dir} \
            --outdir {params.outdir} \
            --locustag {wildcards.sample}_{wildcards.batch} \
            {input.assembly} \
            --force \
            || echo "metaerg.pl exited non-zero (may be output_report.pl failure)"

        GFF_MASTER="{params.outdir}/data/master.gff.txt"
        GFF_ALL="{params.outdir}/data/all.gff"

        if [ -f "$GFF_MASTER" ]; then
            echo "MetaERG completed successfully (master.gff.txt)"
            cp "$GFF_MASTER" {output.gff}
        elif [ -f "$GFF_ALL" ]; then
            echo "MetaERG: using all.gff (output_report.pl known failure)"
            cp "$GFF_ALL" {output.gff}
        else
            echo "ERROR: neither master.gff.txt nor all.gff found under {params.outdir}/data/" >&2
            exit 1
        fi

        # FFN: prefer data/, fall back to tmp/
        for loc in "{params.outdir}/data/cds.ffn" "{params.outdir}/tmp/cds.ffn"; do
            if [ -f "$loc" ]; then
                cp "$loc" {output.ffn}
                break
            fi
        done
        if [ ! -f "{output.ffn}" ]; then
            echo "ERROR: cds.ffn not found in data/ or tmp/" >&2; exit 1
        fi

        # FAA: prefer data/, fall back to tmp/
        for loc in "{params.outdir}/data/cds.faa" "{params.outdir}/tmp/cds.faa"; do
            if [ -f "$loc" ]; then
                cp "$loc" {output.faa}
                break
            fi
        done
        if [ ! -f "{output.faa}" ]; then
            echo "ERROR: cds.faa not found in data/ or tmp/" >&2; exit 1
        fi
        rm -f {params.outdir}/stdin.part_*.tar.gz
        echo "Completed annotation for {wildcards.sample}/{wildcards.batch}"
        """


def get_batch_gffs(wildcards):
    return expand(
        "{sample}/annotation/{batch}/data/either_all_or_master.gff",
        sample=wildcards.sample,
        batch=get_batches(wildcards),
    )


def get_batch_ffns(wildcards):
    return expand(
        "{sample}/annotation/{batch}/data/cds.ffn",
        sample=wildcards.sample,
        batch=get_batches(wildcards),
    )


def get_batch_faas(wildcards):
    return expand(
        "{sample}/annotation/{batch}/data/cds.faa",
        sample=wildcards.sample,
        batch=get_batches(wildcards),
    )


rule join_metaerg_outputs:
    input:
        gff=get_batch_gffs,
        ffn=get_batch_ffns,
    output:
        gff="{sample}/annotation/{sample}_metaerg.gff",  
        ffn="{sample}/annotation/{sample}_metaerg.ffn", 
    container:
        config["docker_seqkit"]
    shell:
        """
        # Concatenate GFF files, keeping only the header from the first file
        grep "^##" {input.gff[0]} > {output.gff}
            for f in {input.gff}; do
        grep -v "^##" "$f" >> {output.gff}
            done

        # Concatenate FFN files (no header stripping needed for FASTA)
        for f in {input.ffn}; do
            cat "$f" >> {output.ffn}
        done
        """


rule antismash:
    container:
        config["docker_antismash"]
    input:
        assembly=lambda wc: config["assembly"][wc.sample],
        gff="{sample}/annotation/{sample}_metaerg.gff",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16 * 1024,
        runtime=6 * 60,
    threads: 16
    log:
        o="{sample}/logs/antismash.log",
    output:
        gbk="{sample}/antismash/{sample}_antismash.gbk",
        outdir=directory("{sample}/antismash/antismash_results/"),
    shell:
        # --genefinding-tool none: ignore contigs without genes
        # https://www.biostars.org/p/9539337/
        # Note: antismash can legitimately fail on very fragmented assemblies;
        # we touch the output so downstream rules are not blocked.
        """
        set +e -x
        antismash \
            --cpus {threads} \
            --allow-long-headers \
            --output-dir {wildcards.sample}/antismash/antismash_results/ \
            {input.assembly} \
            --genefinding-gff {input.gff} \
            --verbose \
            --logfile {log.o} \
            --genefinding-tool none
        exitcode=$?
        if [ $exitcode -ne 0 ]; then
            echo "antismash exited $exitcode; assembly may be too fragmented"
            touch {output.gbk}
        else
            cat {wildcards.sample}/antismash/antismash_results/*.gbk > {output.gbk}
        fi
        """


rule tabulate_antismash:
    input:
        gbk="{sample}/antismash/{sample}_antismash.gbk", 
    output:
        tab="{sample}/antismash/{sample}_antismash.tab",
    container:
        config["docker_biopython"]
    script:
        "../scripts/parse_antismash_gbk.py"


rule annotate_abricate:
    input:
        assembly=lambda wc: config["assembly"][wc.sample],
    output:
        out="{sample}/abricate/{tool}.tab",
    container:
        config["docker_abricate"]
    resources:
        mem_mb=4000,
    shell:
        """
        abricate --db {wildcards.tool} {input.assembly} > {output.out}
        """


rule annotate_AMR:
    input:
        assembly=lambda wc: config["assembly"][wc.sample],
    output:
        amr="{sample}/amrfinder/{sample}_amrfinder.tab",
    resources:
        mem_mb=4000,
        runtime=3 * 60,
    container:
        config["docker_amrfinder"]
    shell:
        """
        amrfinder -n {input.assembly} --plus > {output.amr}
        """


def get_annotate_cazi_runtime(wildcards, attempt):
    return attempt * 3.5 * 60


def get_annotate_cazi_memory(wildcards, attempt):
    return attempt * 4 * 1024


rule annotate_CAZI_split:
    input:
        faa="{sample}/annotation/{batch}/data/cds.faa", 
        gff="{sample}/annotation/{batch}/data/either_all_or_master.gff",
    output:
        overview="{sample}/cazi/batches/{batch}/overview.txt", 
        substrate="{sample}/cazi/batches/{batch}/substrate.out",
        cgc="{sample}/cazi/batches/{batch}/cgc.out",
    params:
        cazi_db=config["cazi_db"],
        outdir=lambda wc: f"{wc.sample}/cazi/batches/{wc.batch}",
    resources:
        mem_mb=get_annotate_cazi_memory,
        runtime=get_annotate_cazi_runtime,
    container:
        config["docker_dbcan"]
    threads: 2
    shell:
        # Pre-touch outputs so Snakemake is satisfied even when run_dbcan
        # produces no hits (empty HMMER results on sparse batches).
        # set -e is active for setup; only run_dbcan itself is wrapped with ||.
        """
        set -euo pipefail
        touch {output.overview} {output.substrate} {output.cgc}
        run_dbcan {input.faa} protein \
            --out_dir {params.outdir}/ \
            -t all \
            --db_dir /app/db \
            -c {input.gff} \
            --cgc_substrate \
            --dia_cpu {threads} \
            --hmm_cpu {threads} \
            --tf_cpu {threads} \
            --stp_cpu {threads} \
            --dbcan_thread {threads} \
            || echo "WARNING: run_dbcan produced no output for {wildcards.sample}/{wildcards.batch}"
        """


rule build_gene_index:
    """Build bowtie2 index from merged gene sequences."""
    input:
        genes="{sample}/annotation/{sample}_metaerg.ffn",
    output:
        multiext(
            "{sample}/annotation/{sample}_genes_index",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
            ".rev.1.bt2", ".rev.2.bt2",
        ),
    threads: 8
    container:
        config["docker_bowtie2"]
    shell:
        """
        bowtie2-build --threads {threads} {input.genes} {wildcards.sample}/annotation/{wildcards.sample}_genes_index
        """


rule align_reads_to_all_genes:
    """Align paired-end reads to the merged gene catalogue."""
    input:
        idx=multiext(
            "{sample}/annotation/{sample}_genes_index",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
            ".rev.1.bt2", ".rev.2.bt2",
        ),
        r1=lambda wc: config["R1"][wc.sample],
        r2=lambda wc: config["R2"][wc.sample],
    output:
        bam=temp("{sample}/annotation/{sample}_aligned_reads.bam"),
        bai=temp("{sample}/annotation/{sample}_aligned_reads.bam.bai"),
    threads: 16
    resources:
        mem_mb=16 * 1024,
        runtime=120,
    container:
        config["docker_bowtie2"]
    shell:
        """
        bowtie2 --threads {threads} \
            -1 {input.r1} -2 {input.r2} \
            -x {wildcards.sample}/annotation/{wildcards.sample}_genes_index \
            | samtools view -@ {threads} -Sb \
            | samtools sort -o {output.bam} -@ {threads}
        samtools index {output.bam}
        """


rule calculate_gene_coverage:
    """Per-base 5'-end read depth across all genes (bedtools genomecov)."""
    input:
        bam="{sample}/annotation/{sample}_aligned_reads.bam",
        bai="{sample}/annotation/{sample}_aligned_reads.bam.bai",
    output:
        coverage=temp("{sample}/annotation/{sample}_gene_coverage.txt"),
    container:
        config["docker_bedtools"]
    shell:
        """
        bedtools genomecov -d -5 -ibam {input.bam} > {output.coverage}
        """


def get_cazi_overviews(wildcards):
    return expand(
        "{sample}/cazi/batches/{batch}/overview.txt",
        sample=wildcards.sample,
        batch=get_batches(wildcards),
    )


def get_cazi_substrates(wildcards):
    return expand(
        "{sample}/cazi/batches/{batch}/substrate.out",
        sample=wildcards.sample,
        batch=get_batches(wildcards),
    )


def get_cazi_cgcs(wildcards):
    return expand(
        "{sample}/cazi/batches/{batch}/cgc.out",
        sample=wildcards.sample,
        batch=get_batches(wildcards),
    )

rule flagstat:
    input:
        bam="{sample}/annotation/{sample}_aligned_reads.bam",
    output:
        flagstat="{sample}/annotation/{sample}_flagstat.txt",
    container:
        config["docker_bowtie2"]
    shell:
        "samtools flagstat {input.bam} > {output.flagstat}"

rule join_CAZI:
    """Merge per-batch CAZI results and compute RPM against global coverage."""
    input:
        overview=get_cazi_overviews,
        substrate=get_cazi_substrates,
        cgc=get_cazi_cgcs,
        coverage="{sample}/annotation/{sample}_gene_coverage.txt",
        flagstat="{sample}/annotation/{sample}_flagstat.txt",
    output:
        overview="{sample}/cazi/{sample}_cazi_overview.txt",
        substrate="{sample}/cazi/{sample}_cazi_substrate.out",
        cgc="{sample}/cazi/{sample}_cazi_cgc.out", 
        rpm="{sample}/cazi/{sample}_annotated_cazymes_RPM.tsv",
    conda:
        "../envs/annotate_output_parse.yaml"
    script:
        "../scripts/join_and_generate_RPM.py"
