import os
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


configfile: os.path.join(str(workflow.basedir), "../../config/config.yaml")


envvars:
    "TMPDIR",

SHARDS = make_shard_names(config["nshards"])


localrules:
    all,


quast_outputs = expand(
    ["quast/quast_{sample}/report.pdf",
    "quast/quast_{sample}/transposed_report.tsv"],
    sample=config["sample"]
)
all_inputs = quast_outputs


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")
    if config["assembler"].lower() == "spades":
        if len(config["R1"]) > 1:
            print("WARNING: Concatenating multiple inputs into a single paired library as metaspades does not support multiple libraries")
        print("Running SPAdes")
    else:
        print("Running megahit")


if config["assembler"].lower() == "spades":
    assemblies = expand(
        ["spades_{sample}.assembly.fasta",
        "spades_{sample}_metaviral/scaffolds.fasta"],
        sample=config["sample"]
    )
    assemblies_labels = ",".join(
        [f"metaspades_{s},metaviralspades_{s}" for s in config["sample"]]
    )
    all_inputs.extend(expand("{sample}.cleaned_assembly_files",
    sample=config["sample"]))
    all_inputs.extend(assemblies)
# elif config["assembler"].lower() == "both": # should this be enabled?
#     assemblies = [f"spades_{config['sample']}.assembly.fasta", f"spades_{config['sample']}_metaviral/scaffolds.fasta", f"megahit_{config['sample']}.assembly.fasta"]
#     assemblies_labels = ",".join([f"metaspades_{config['sample']}",f"metaviralspades_{config['sample']}",f"megahit_{config['sample']}"])
#     all_inputs.append(f"{config['sample']}.cleaned_assembly_files")
#     all_inputs.extend(assemblies)
else:
    assemblies = expand("megahit_{sample}.assembly.fasta", sample=config["sample"])
    assemblies_labels = ",".join([f"megahit_{s}" for s in config["sample"]])
    all_inputs.extend(assemblies)


rule all:
    input:
        all_inputs,


#
if len(config["sample"]) == 1 and len(config["R1"][config["sample"][0]]) == 1:
    input_R1 = config["R1"][config["sample"][0]]
    input_R2 = config["R2"][config["sample"][0]]
else:
    input_R1 = expand("concatenated/{sample}_R1.fastq.gz", sample=config["sample"])
    input_R2 = expand("concatenated/{sample}_R2.fastq.gz", sample=config["sample"])


# Utils Module
module utils:
    snakefile: 
       "utils.smk"
    config: 
       config
    skip_validation: 
       True


'''
#Jul/16/2025: I commented this because this function is defined in many rules. It was given a trouble in multisample preprocess.
#I can do `use rule * from assembly exclude_rules: concat_lanes_fix_names` in Snakefile to try to fix it.
use rule concat_lanes_fix_names from utils as utils_concat_lanes_fix_names with:
    input:
        fq=get_concat_input_multisample,
    output:
        fq=temp("concatenated/{sample}_R{rd}.fastq.gz"),
    log:
        e="logs/concat_lanes_fix_names_{sample}_R{rd}.e",
'''



rule megahit:
    input:
        R1=lambda wc: config["R1"][wc.sample],
        R2=lambda wc: config.get("R2", {}).get(wc.sample, [])
    output:
        outdir=directory("megahit_{sample}"),
        assembly="megahit_{sample}.assembly.fasta"
    container:
        config["docker_megahit"]
    resources:
        mem_mb=64000,
        runtime=24 * 60
    threads: 64
    params:
        input_string=lambda wc, input: (
            "-1 " + ",".join(input.R1) + " -2 " + ",".join(input.R2)
            if input.R2 else "-r " + ",".join(input.R1)
        )
    shell:
        """
        mkdir -p ${{TMPDIR}}/megahit_{wildcards.sample}/
        megahit {params.input_string} \
            --out-dir megahit_{wildcards.sample}/ \
            --out-prefix {wildcards.sample} \
            --tmp-dir ${{TMPDIR}}/megahit_{wildcards.sample}/ \
            --memory $((64000 * 1024)) \
            --num-cpu-threads {threads}
        rm -r ${{TMPDIR}}/megahit_{wildcards.sample}/
        mv megahit_{wildcards.sample}/{wildcards.sample}.contigs.fa {output.assembly}
        """


rule SPAdes_run:
    # TODO: add in params for read length to experiment with larger kmers than default
    input:
        get_config_inputs_multisample
    output:
        assembly="spades_{sample}.assembly.fasta",
        graph="spades_{sample}.assembly_graph.gfa",
        spades_log="spades_{sample}/spades.log",
    container:
        config["docker_spades"]
    conda:
        "../envs/spades.yaml"
    resources:
        mem_mb=64000,
        #mem_mb=lambda wc, attempt, input: attempt
        #* (max(sum(Path(f).stat().st_size for f in input.R1 + input.R2) // 1000000, 1024) * 20),
        runtime=48 * 60,
    params:
        input_string=lambda wc, input: "-1 " + ",".join(input.R1) + " -2 " + ",".join(input.R2),
    threads: 64
    log:
        e="logs/spades_{sample}.log",
    shell:
        """
        spades.py \
            -1 {input.R1} \
            -2 {input.R2} \
            -t {threads} \
            --meta \
            -o spades_{wildcards.sample} \
            --tmp-dir $TMPDIR \
            -m $(({resources.mem_mb}/1024)) \
            2> {log.e}
        mv spades_{wildcards.sample}/scaffolds.fasta {output.assembly}
        mv spades_{wildcards.sample}/assembly_graph_with_scaffolds.gfa {output.graph}
        """


rule viral_SPAdes_run:
    # max_kmer https://github.com/ablab/spades/discussions/1188
    # --onlyassembler seems to be neccessary when using assembly graph input
    input:
        get_config_inputs_multisample,
        #R1=lambda wc: config["R1"][wc.sample],
        #R2=lambda wc: config["R2"][wc.sample],
        assembly_graph="spades_{sample}.assembly_graph.gfa",
    output:
        assembly="spades_{sample}_metaviral/scaffolds.fasta",
        spades_log="spades_{sample}_metaviral/spades.log",
    container:
        config["docker_spades"]
    params:
        max_kmer=55,
    conda:
        "../envs/spades.yaml"
    resources:
        mem_mb=lambda wc, attempt, input: attempt
        * (max(sum(Path(f).stat().st_size for f in input.R1 + input.R2) // 1000000, 1024) * 20),
        runtime=48 * 60,
    threads: 64
    log:
        e="logs/spades_{sample}_metaviral.log",
    shell:
        """
        spades.py \
            --metaviral \
            -1 {input.R1} \
            -2 {input.R2} \
            --assembly-graph {input.assembly_graph} \
            --only-assembler \
            -t {threads} \
            -o spades_{wildcards.sample}_metaviral/ \
            -m $(({resources.mem_mb}/1024)) \
            -k {params.max_kmer}  \
            2> {log.e}
        # deal with missing scaffolds file if no viruses recovered
        if grep -q "No complete extrachromosomal contigs assembled" "spades_{wildcards.sample}_metaviral/spades.log"; then
            touch {output.assembly}
        fi
        """


rule quast_run:
    """ see http://quast.sourceforge.net/docs/manual.html
    NOTE: giving --blast-db {input.blast_16s_db_nsq} as the NCBI 16s db
    results in silent errors because their automatic genome retrieval
    uses names not accessions :(
    so instead we allow it to use its own SILVA one
    I tried downloading their silva and preprocessing it but
    ran into additional errors, presumably due to the OLD version
    of BLAST used by quast...
    TODO: make --reference work with chocophlan db somehow?

    As of 2022-11-09 we switched to non-meta quast for runtime issues
    with pulling genomes
    """
    input:
        assembly=lambda wc: ([f"spades_{wc.sample}.assembly.fasta", f"spades_{wc.sample}_metaviral/scaffolds.fasta"]
                if config["assembler"].lower() == "spades"
                else [f"megahit_{wc.sample}.assembly.fasta"]),
        blast_16s_db_nsq=config["blast_16s_db_nsq"],
    output:
        report_tsv="quast/quast_{sample}/transposed_report.tsv",
        report="quast/quast_{sample}/report.pdf",
    params:
        dir="quast/quast_{sample}/",
        labels=lambda wc: (
                f"metaspades_{wc.sample},metaviralspades_{wc.sample}"
                if config["assembler"].lower() == "spades"
                else f"megahit_{wc.sample}"),
    threads: 16
    container:
        config["docker_quast"]
    conda:
        "../envs/spades.yaml"
    log:
        e="logs/quast_{sample}.e",
        o="logs/quast_{sample}.o",
    resources:
        mem_mb=8 * 1024,
        runtime=3 * 60,
    shell:
        """
        quast.py \
          -o {params.dir} \
          --split-scaffolds \
          --threads {threads} \
          {input.assembly} \
          --ambiguity-usage all \
          --min-contig 100 \
          --no-snps \
          --no-icarus \
          --labels {params.labels} \
          > {log.o} 2> {log.e}
        """


rule clean_up:
    """"{sample}_metaerg.gff" is used as an input to ensure
    that step is done before we clean.
    """
    input:
        agg_files=lambda wc: [f"quast/quast_{wc.sample}/report.pdf"],
        spades_log="spades_{sample}/spades.log",
    output:
        touch("{sample}.cleaned_assembly_files"),
    shell:
        """
        # We dont need any of the spades intermediate files (corrected reads, temp, per-kmer assembly)
        ls spades_{wildcards.sample}/
        find spades_{wildcards.sample}/  -type f | xargs --no-run-if-empty rm
        """

