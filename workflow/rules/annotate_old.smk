import os
import sys
from shutil import rmtree
import glob
from math import ceil


include: "common.smk"


configfile: os.path.join(workflow.basedir, "../../config/config.yaml")

wildcard_constraints:
    sample="[^/]+",  # sample cannot contain forward slashes
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

SAMPLES = config["sample"] if isinstance(config["sample"], list) else [config["sample"]]



abricates = expand(
    "{sample}_abricate_{tool}.tab",
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

outputs = (expand("{sample}_antismash.gbk", sample=config['sample']) +
    expand("{sample}_antismash.tab", sample=config['sample']) +
    expand("{sample}_amrfinder.tab", sample=config['sample']) +
    expand("{sample}_cazi_overview.txt", sample=config['sample']) +
    expand("{sample}_cazi_substrate.out", sample=config['sample']) +
    expand("{sample}_annotated_cazymes_RPM.tsv", sample=config['sample']) +
    abricates)
if config.get("check_contigs", False):
outputs += expand("{sample}_metaerg.gff", sample=config['sample'])

rule all:
    input:
        outputs,


checkpoint split_assembly:
    """Split assembly into chunks for parallel processing"""
    input:
        assembly=lambda wc: config["assembly"][wc.sample] if isinstance(config["assembly"], dict) else config["assembly"],
    output:
        splitdir=directory("tmp/{sample}"),
        #assembly="tmp/{sample}/assembly.fasta",
        done="tmp/{sample}/.split_done",
    params:
        outdir=lambda wc: os.path.abspath("tmp/" + wc.sample),
        #assembly_copy=lambda wc: os.path.abspath("tmp/" + wc.sample + "/assembly.fasta"),
        minlen=config["contig_annotation_thresh"],
        nseqs=config.get("chunk_size_contigs", 200),
        nbatches=lambda wc, input: calculate_nbatches(
            str(input.assembly),
            config["contig_annotation_thresh"],
            config.get("chunk_size_contigs", 200)
        ),
    container:
        "docker://pegi3s/seqkit:2.3.0"
    threads: 4
    resources:
        mem_mb=8000,
    log:
        e="logs/split_assembly_{sample}.e",
        o="logs/split_assembly_{sample}.o",
    shell:
        """
        set -e
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        seqkit shuffle {input.assembly} --two-pass \
        | seqkit seq --min-len {params.minlen} --threads {threads} \
        | seqkit split2 --by-part {params.nbatches} --out-dir {params.outdir} --force > "{log.o}" 2>> "{log.e}"
       
        touch {output.done}
        """


def calculate_nbatches(assembly_path, minlen, nseqs):
    """Calculate number of batches needed for splitting assembly"""
    from math import ceil
    
    n_filtered = count_contigs_with_minlen(assembly_path, minlen)
    if n_filtered <= 0:
        nbatches = 1
    else:
        nbatches = min(ceil(n_filtered / nseqs), n_filtered)
    
    print(f"Detected {n_filtered} contigs >= {minlen} bp; splitting into {nbatches} parts")
    return nbatches

rule annotate_orfs:
    container:
        config["docker_metaerg"]
    input:
        #splitdir="tmp/{sample}",
        assembly="tmp/{sample}/{batch}.fasta",
    output:
        gff="annotation/{sample}/annotation_{batch}/data/either_all_or_master.gff",
        ffn="annotation/{sample}/annotation_{batch}/data/cds.ffn",
        faa="annotation/{sample}/annotation_{batch}/data/cds.faa",
    resources:
        mem_mb=8 * 1024,
        runtime=45,
    threads: 4
    params:
        metaerg_db_dir=config["metaerg_db_dir"],
    shell:
        """
        # turn off strict so we don't fail even if we have the gff file.
        # currently the output_report.pl script is failing
        # see issues https://github.com/xiaoli-dong/metaerg/pull/38 and
        # https://github.com/xiaoli-dong/metaerg/issues/12
        set +e
        metaerg.pl --cpus {threads} --dbdir {params.metaerg_db_dir} --outdir annotation/{wildcards.sample}/annotation_{wildcards.batch} --locustag {wildcards.sample}_{wildcards.batch} {input.assembly} --force || echo "Finished running Metaerg"
        # if metaerg successfully packaged everything up
        if [ -f "annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/master.gff.txt" ]
        then
            echo "MetaERG completed successfully"
            cp annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/master.gff.txt {output.gff}
        elif [ -f "annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/all.gff" ]; then
            echo "MetaERG annotation succeeded but output_report.pl failed (known issue)"
            cp annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/all.gff {output.gff}
        else
            # if it successed but failed at output_report.pl, no need to do anything
            echo "sample likely failed at output_report.pl but gff should be present"
            mv annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/all.gff {output.gff}
        fi

        # Check for FFN file in data directory first, then tmp directory
        if [ -f "annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/cds.ffn" ]; then
            cp annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/cds.ffn {output.ffn}
        elif [ -f "annotation/{wildcards.sample}/annotation_{wildcards.batch}/tmp/cds.ffn" ]; then
            echo "Copying cds.ffn from tmp directory (output_report.pl failed)"
            cp annotation/{wildcards.sample}/annotation_{wildcards.batch}/tmp/cds.ffn {output.ffn}
        else
            echo "ERROR: cds.ffn not found in data or tmp directory"
            exit 1
        fi
        
        # Check for FAA file in data directory first, then tmp directory
        if [ -f "annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/cds.faa" ]; then
            cp annotation/{wildcards.sample}/annotation_{wildcards.batch}/data/cds.faa {output.faa}
        elif [ -f "annotation/{wildcards.sample}/annotation_{wildcards.batch}/tmp/cds.faa" ]; then
            echo "Copying cds.faa from tmp directory (output_report.pl failed)"
            cp annotation/{wildcards.sample}/annotation_{wildcards.batch}/tmp/cds.faa {output.faa}
        else
            echo "ERROR: cds.faa not found in data or tmp directory"
            exit 1
        fi
        
        echo "Successfully completed annotation for {wildcards.sample}_{wildcards.batch}"
        """

def get_batch_gffs(wildcards):
    """Aggregate all GFF files for a sample after checkpoint completes"""
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    batch_names = [os.path.basename(b).replace('.fasta', '') for b in batches]
    
    return expand(
        "annotation/{sample}/annotation_{batch}/data/either_all_or_master.gff",
        sample=wildcards.sample,
        batch=batch_names
    )


def get_batch_ffns(wildcards):
    """Aggregate all FFN files for a sample after checkpoint completes"""
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    batch_names = [os.path.basename(b).replace('.fasta', '') for b in batches]

    return expand(
        "annotation/{sample}/annotation_{batch}/data/cds.ffn",
        sample=wildcards.sample,
        batch=batch_names
    )


rule join_metaerg_outputs:
    input:
        gff=get_batch_gffs,
        ffn=get_batch_ffns,
    output:
        gff="{sample}_metaerg.gff",
        ffn="{sample}_metaerg.ffn",
    container:
        config["docker_seqkit"]
    shell:
        """
        # deal with header
        head -n 1 {input.gff[0]} > {output.gff}
        for f in {input.gff}
        do
            tail -n+2 $f >> {output.gff}
        done
        for f in {input.ffn}
        do
            cat $f >> {output.ffn}
        done
        """

rule antismash:
    # note that we don't require the web index as antismash can fail on small samples.
    container:
        config["docker_antismash"]
    input:
        assembly=lambda wc: config["assembly"][wc.sample],
        #assembly=config["assembly"],
        gff="{sample}_metaerg.gff",
    resources:
        mem_mb=16 * 1024,
        runtime=6 * 60,
    threads: 16
    log:
        o="logs/antismash_{sample}.log",
    output:
        gbk="{sample}_antismash.gbk",
        outdir=directory("antismash_{sample}"),
    shell:
        """
        set +e -x
        antismash --cpus {threads} --allow-long-headers \
            --output-dir antismash_{wildcards.sample} {input.assembly} \
            --genefinding-gff {input.gff} --verbose --logfile {log.o} \
            --genefinding-tool none  # ignore contigs without genes: https://www.biostars.org/p/9539337/
        exitcode=$?
        if [ ! $exitcode -eq 0 ]
        then
            echo "Error running antismash; this can occur if the assembly is very fragmented/poor"
            touch {output.gbk}

        else
            cat antismash_{wildcards.sample}/*.gbk > {output.gbk}
        fi
        """


rule tabulate_antismash:
    input:
        gbk="{sample}_antismash.gbk",
    output:
        tab="{sample}_antismash.tab",
    container:
        config["docker_biopython"]
    script:
        "../scripts/parse_antismash_gbk.py"


rule annotate_abricate:
    input:
        #assembly=config["assembly"]
        assembly=lambda wc: config["assembly"][wc.sample],
    output:
        out="{sample}_abricate_{tool}.tab",
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
        #assembly=config["assembly"]
        assembly=lambda wc: config["assembly"][wc.sample],
    output:
        amr="{sample}_amrfinder.tab",
    resources:
        mem_mb=4000,
        runtime=3 * 60,
    container:
        config["docker_amrfinder"]
    shell:
        """
        amrfinder -n {input.assembly}  --plus  > {output.amr}
        """

def get_annotate_cazi_runtime(wildcards, attempt):
    return attempt * 3.5 * 60


def get_annotate_cazi_memory(wildcards, attempt):
    return attempt * 4 * 1024


def get_batch_faas(wildcards):
    """Get all FAA files for CAZI annotation"""
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    batch_names = [os.path.basename(b).replace('.fasta', '') for b in batches]
    return batch_names



rule annotate_CAZI_split:
    input:
        faa="annotation/{sample}/annotation_{batch}/data/cds.faa",
        gff="annotation/{sample}/annotation_{batch}/data/either_all_or_master.gff",
    output:
        overview="cazi_db_scan/{sample}/{batch}/overview.txt",
        substrate="cazi_db_scan/{sample}/{batch}/substrate.out",
        cgc="cazi_db_scan/{sample}/{batch}/cgc.out",
    params:
        cazi_db=config["cazi_db"],
        contig_annotation_thresh=config["contig_annotation_thresh"],
    resources:
        mem_mb=get_annotate_cazi_memory,
        runtime=get_annotate_cazi_runtime,
    container:
        config["docker_dbcan"]
    threads: 2
    shell:
        """
        # turn off strict so we don't fail on empty hmmer outputs within run_dbcan
        set -e
        # touch this file so it exists even if this fails
        # if you can figure out checkpoints that can deal with missing files,
        # you win!
        touch {output.overview}
        touch {output.substrate}
        touch {output.cgc}
        run_dbcan {input.faa} protein --out_dir cazi_db_scan/{wildcards.sample}/{wildcards.batch}/ -t all --db_dir /app/db -c {input.gff} --cgc_substrate --dia_cpu {threads} --hmm_cpu {threads} --tf_cpu {threads} --stp_cpu {threads} --dbcan_thread {threads} ||  echo "WARNING: no output from this split!"
        """

rule align_annotated_genes:
    input:
        ffn="annotation/{sample}/annotation_{batch}/data/cds.ffn",
        r1=lambda wc: config["R1"][wc.sample],
        r2=lambda wc: config["R2"][wc.sample],
    output:
        bamfile="annotation/{sample}/annotation_{batch}/aligned_reads.bam",
    container:
        config["docker_bowtie2"]
    resources:
        mem_mb=16 * 1024,
        runtime=get_annotate_cazi_runtime,
        threads=16,
        cores=16,
    params:
        bowtie_dir="annotation/{sample}/annotation_{batch}/bowtie",
        bowtie_index="annotation/{sample}/annotation_{batch}/bowtie/bowtie2_index",
    shell:
        """
        mkdir -p {params.bowtie_dir}
        bowtie2-build \
            --threads {resources.threads} \
            {input.ffn} \
            {params.bowtie_index}
        bowtie2 --threads {resources.threads} -1 {input.r1} -2 {input.r2} -x {params.bowtie_index}  | samtools view -@ {resources.threads} -Sb | samtools sort -o {output.bamfile} -@ {resources.threads}
        """


rule seqkit_annotate_ffn:
    input:
        ffn="annotation/{sample}/annotation_{batch}/data/cds.ffn",
    output:
        length_file="annotation/{sample}/annotation_{batch}/seqkit.length",
        bed_file="annotation/{sample}/annotation_{batch}/seqkit.bed",
    container:
        config["docker_seqkit"]
    shell:
        """
        seqkit fx2tab -l -n -i {input.ffn} | awk '{{print $1"\t"$2}}' > {output.length_file}
        seqkit fx2tab -l -n -i {input.ffn} | awk '{{print $1"\t"0"\t"$2}}' > {output.bed_file}
        """


rule bedtools_coverage:
    input:
        length_file="annotation/{sample}/annotation_{batch}/seqkit.length",
        bed_file="annotation/{sample}/annotation_{batch}/seqkit.bed",
        bamfile="annotation/{sample}/annotation_{batch}/aligned_reads.bam",
    output:
        coverage="annotation/{sample}/annotation_{batch}/annotated_gene_coverage.txt",
    container:
        config["docker_bedtools"]
    shell:
        """
        bedtools genomecov -5 -ibam {input.bamfile} > {output.coverage}
        """

rule create_RPM_counts:
    input:
        coverage="annotation/{sample}/annotation_{batch}/annotated_gene_coverage.txt",
        overview="cazi_db_scan/{sample}/{batch}/overview.txt",
        substrate="cazi_db_scan/{sample}/{batch}/substrate.out",
        cgc="cazi_db_scan/{sample}/{batch}/cgc.out",
        r1=lambda wc: config["R1"][wc.sample]
    output:
        rpm_file="cazi_db_scan/{sample}/{batch}/annoted_cazymes_RPM.tsv",
    conda:
        "../envs/annotate_output_parse.yaml"
    script:
        "../scripts/generate_RPM_annotation_files.py"



def get_cazi_overviews(wildcards):
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    batch_names = [os.path.basename(b).replace('.fasta', '') for b in batches]
    return expand("cazi_db_scan/{sample}/{batch}/overview.txt", sample=wildcards.sample, batch=batch_names)


def get_cazi_substrates(wildcards):
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    batch_names = [os.path.basename(b).replace('.fasta', '') for b in batches]
    return expand("cazi_db_scan/{sample}/{batch}/substrate.out", sample=wildcards.sample, batch=batch_names)


def get_cazi_cgcs(wildcards):
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    batch_names = [os.path.basename(b).replace('.fasta', '') for b in batches]
    return expand("cazi_db_scan/{sample}/{batch}/cgc.out", sample=wildcards.sample, batch=batch_names)


def get_cazi_rpms(wildcards):
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output.splitdir
    batches = glob.glob(f"{checkpoint_output}/stdin.part_*.fasta")
    batch_names = [os.path.basename(b).replace('.fasta', '') for b in batches]
    return expand("cazi_db_scan/{sample}/{batch}/annoted_cazymes_RPM.tsv", sample=wildcards.sample, batch=batch_names)

rule join_CAZI:
    input:
        overview=get_cazi_overviews,
        substrate=get_cazi_substrates,
        cgc=get_cazi_cgcs,
        rpm=get_cazi_rpms,
    output:
        overview="{sample}_cazi_overview.txt",
        substrate="{sample}_cazi_substrate.out",
        cgc="{sample}_cazi_cgc.out",
        rpm="{sample}_annotated_cazymes_RPM.tsv",
    shell:
        """
        join_files(){{
            output_file=$1
            echo $output_file
            # deal with header
            head -n 1 $2 > $output_file
            # This weird for loop handles the fact we pass the infiles as
            # a list - (which doesn't work well in bash):
            for (( i=2; i <= "$#"; i++ )); do
                tail -n+2 ${{!i}} >> $output_file
            done
        }}
        
        join_files {output.overview} {input.overview}
        join_files {output.substrate} {input.substrate}
        join_files {output.cgc} {input.cgc}
        join_files {output.rpm} {input.rpm}
        """

