import os
import json
import yaml
import shutil

from pathlib import Path


include: "common.smk"


# the str is needed as when running from github workflow.current_basedir is a Githubfile, not a string or a path so os.path objects
configfile: os.path.join(str(workflow.current_basedir), "../../config/config.yaml")


envvars:
    "TMPDIR",


SHARDS = make_shard_names(config["nshards"])


onstart:
    with open("config_used.yaml", "w") as outfile:
        yaml.dump(config, outfile)

    if not os.path.exists("logs"):
        os.makedirs("logs")


localrules:
    all,

BINNING_TOOLS = ["concoct", "metabat2", "maxbin2"]

def sample_dir(sample):
    return f"{sample}"

def refined_dir(sample):
    return (
        f"{sample_dir(sample)}/refined_binning/"
        f"metawrap_{config['metawrap_compl_thresh']}_"
        f"{config['metawrap_contam_thresh']}_bins"
    )

def mag_dir(sample):
    return refined_dir(sample)

#binstats_all = expand(
#    "metawrap/rawbinning_{sample}/{tool}/{tool}_bins/{tool}.done",
#    sample=config["sample"],
#    tool=BINNING_TOOLS,)

refined_binstats_each = expand(
    "{sample}/refined_binning/{tool}_bins.stats",
    sample=config["sample"],
    tool=BINNING_TOOLS,)

refined_stats = expand(
    f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats',
    sample=config["sample"],)

stats_mqc = expand(
    f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats_mqc.tsv',
    sample=config["sample"],)

contigs = expand(
    f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.contigs',
    sample=config["sample"],)

covermreports = expand(
    f'{{sample}}/coverm/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.coverage_mqc.tsv',
    sample=config["sample"],)
gtdbk = expand("{sample}/gtdbtk/", sample=config["sample"])
prokka = expand("{sample}/prokka.done", sample=config["sample"])

rule all:
    input:
        refined_stats,
        refined_binstats_each,
        stats_mqc,
        contigs,
        covermreports,
        gtdbk,
        prokka,
        expand("{sample}/logs/cleanup_rawbinning.done", sample=config["sample"]),


rule unzip_rename_fastq_for_metawrap:
    input:
        R1=lambda wc: config["R1"][wc.sample],
        R2=lambda wc: config["R2"][wc.sample],
    output:
        R1=temp("tmp_{sample}_1.fastq"),
        R2=temp("tmp_{sample}_2.fastq"),
    shell:
        """
        zcat {input.R1} > {output.R1}
        zcat {input.R2} > {output.R2}
        """


rule metawrap_binning:
    """
    """
    input:
        R1="tmp_{sample}_1.fastq",
        R2="tmp_{sample}_2.fastq",
        assembly=lambda wc: config["assembly"][wc.sample],
    output:
        stats=temp("metawrap/rawbinning_{sample}/{tool}/{tool}_bins/{tool}.done"),
    params:
        outdir=lambda wc, output: os.path.dirname(os.path.dirname(output.stats)),
    container:
        config["docker_metawrap"]
    threads: 64
    resources:
        mem_mb=32 * 1024,
        runtime=12 * 60,
    shell:
        """
        echo "tool is: '{wildcards.tool}'"
        echo "outdir is: '{params.outdir}'"
        metawrap binning -o {params.outdir} -t {threads} -a {input.assembly} --{wildcards.tool} {input.R1} {input.R2}
        touch {output.stats}
        """

def get_binstats(wildcards):
    return expand(
        "metawrap/rawbinning_{sample}/{tool}/{tool}_bins/{tool}.done",
        sample=wildcards.sample,
        tool=BINNING_TOOLS,
    )


checkpoint metawrap_refine_binning:
    """The names for the params binput_dirs is due to the crazy naming of the output of metawrap.
    The only consistant file we can use as a trigger is the <tool>.done file, but the refine module needs the dir beneath it under a path like
    metawrap/rawbinning_473/concoct/concoct_bins/concoct_bins
    """
    input:
        binputs=get_binstats,
    output:
        stats=f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats',
        stats_mqc=f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats_mqc.tsv',
        contigs=f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.contigs',
        bin1=f"{{sample}}/refined_binning/{BINNING_TOOLS[0]}_bins.stats",
        bin2=f"{{sample}}/refined_binning/{BINNING_TOOLS[1]}_bins.stats",
        bin3=f"{{sample}}/refined_binning/{BINNING_TOOLS[2]}_bins.stats",
    params:
        outdir=lambda wc, output: os.path.dirname(output.stats),
        binput_dirs=lambda wc, input: [os.path.dirname(x) for x in input.binputs],
        completeness=config["metawrap_compl_thresh"],
        contamination=config["metawrap_contam_thresh"],
        checkm_db=config["checkm_db"],
    # this is a different container that has checkm installed
    container:
        config["docker_metawrap"]
    threads: 32
    # give it 82 gb memory because checkm estimates pplacer will need 40GB per core, and we want this to run reasonably fast
    resources:
        mem_mb=82 * 1024,
        runtime=12 * 60,
    shell:
        """
        export  CHECKM_DATA_PATH={params.checkm_db}
        metawrap bin_refinement -o {params.outdir} -t {threads} -A {params.binput_dirs[0]} -B {params.binput_dirs[1]} -C {params.binput_dirs[2]} -c {params.completeness} -x {params.contamination}
        echo -e "#id: 'metawrap'\n#plot_type: 'table'\n#section_name: 'Bin Refinement'" > {output.stats}_mqc.tsv && cat {output.stats} >> {output.stats}_mqc.tsv
        """

def get_fastqs(config_entry):
    """Always returns a list regardless of whether config value is a string or list"""
    if isinstance(config_entry, str):
        return [config_entry]
    return config_entry


rule coverm:
    """ This calculates bin coverage.  The bin stats file is used as the
    trigger, but we actually want the directory of the bins

    specifying inputs with --coupled and -1 -2 seem to give equivalent results.
    Note that the output only refers tothe forward file as the sample name.
    minimap2 is supposedly faster and more accurate than BWA, so we use that
    as the mapper.
    we return all the available methods as of the time of writing this rule
    except for coverage_histogram which has to run separately
    ("Cannot specify the coverage_histogram method with any other coverage methods")

    --min-covered-fraction is required to be set to 0 for certain cov metrics
     --genome-fasta-extension is fa since thats how metawrap outputs it
    """
    input:
        R1=lambda wc: config["R1"][wc.sample],
        R2=lambda wc: config["R2"][wc.sample],
        stats=f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats',
    output:
        mqc=f'{{sample}}/coverm/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.coverage_mqc.tsv',
        bams=directory(
            f'{{sample}}/coverm/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bams/'),
    params:
        bindir=lambda wc, input: input.stats.replace(".stats", ""),
        fastq_string=lambda wc: " ".join(
              [r1 + " " + r2
                for r1, r2 in zip(
                   get_fastqs(config["R1"][wc.sample]),
                   get_fastqs(config["R2"][wc.sample]),
                )]),
    container:
        config["docker_coverm"]
    threads: 16
    shell:
        """
        coverm genome --genome-fasta-directory {params.bindir} \
          --coupled {params.fastq_string} \
          --mapper minimap2-sr \
          --methods mean relative_abundance trimmed_mean \
            covered_bases variance length count reads_per_base rpkm tpm \
          --output-file {output.mqc}.tmp --threads {threads} \
          --bam-file-cache-directory {output.bams} \
          --min-covered-fraction 0 \
          --genome-fasta-extension fa
        # add in the Multiqc header info
        echo -e "# plot_type: 'table'\n# section_name: 'Bin Coverage Statistics'" > {output.mqc}
        cat {output.mqc}.tmp >> {output.mqc}
        rm {output.mqc}.tmp
        """



#def get_prokka_bins(wildcards):
#    bin_dir = "metawrap/refined_binning_{}/metawrap_{}_{}_{}/".format(
#        config['sample'],
#        config["metawrap_compl_thresh"],
#        config["metawrap_contam_thresh"],
#        "bins"
#    )
#    bins = glob_wildcards(bin_dir + "{bin}.fa").bin
#    return expand("prokka/{sample}/{bin}/",
#                  sample=config['sample'],
#                  bin=bins)


#def mag_dir(sample):
#    return (
#        f"metawrap/refined_binning_{sample}/"
#        f"metawrap_{config['metawrap_compl_thresh']}_"
#        f"{config['metawrap_contam_thresh']}_bins"
#    )

#def get_mag_bins(wildcards):
#    """Dynamically get all bin .fa files for a sample after refinement."""
#    d = mag_dir(wildcards.sample)
#    return sorted(Path(d).glob("*.fa"))

#def get_prokka_outputs(wildcards):
#    """Get all expected Prokka output dirs for a sample."""
#    d = mag_dir(wildcards.sample)
#    bins = [f.stem for f in sorted(Path(d).glob("*.fa"))]
#    return expand(
#        "prokka/{sample}/{bin}/",
#        sample=wildcards.sample,
#        bin=bins,
#    )
#def get_prokka_outputs(wildcards):
#    stats_file = (
#        f"{wildcards.sample}/refined_binning/"
#        f"metawrap_{config['metawrap_compl_thresh']}_"
#        f"{config['metawrap_contam_thresh']}_bins.stats"
#    )
#    with open(stats_file) as f:
#        next(f)  # skip header
#        bins = [line.split("\t")[0] for line in f if line.strip()]
#    return expand(
#        "{sample}/prokka/{bin}/",
#        sample=wildcards.sample,
#        bin=bins,
#    )
def get_prokka_outputs(wildcards):
    # CHANGE: use checkpoint to wait for refinement before reading stats
    checkpoints.metawrap_refine_binning.get(sample=wildcards.sample)
    stats_file = (
        f"{wildcards.sample}/refined_binning/"
        f"metawrap_{config['metawrap_compl_thresh']}_"
        f"{config['metawrap_contam_thresh']}_bins.stats"
    )
    with open(stats_file) as f:
        next(f)
        bins = [line.split("\t")[0] for line in f if line.strip()]
    return expand(
        "{sample}/prokka/{bin}/",
        sample=wildcards.sample,
        bin=bins,
    )
#rule gtdbtk_classify_wf:
#    """
#    Taxonomic classification of MAGs using GTDB-Tk.
#    Input: refined MAG bins from MetaWRAP
#    Output: taxonomic classification results per sample
#    """
#    input:
#        bins = expand("metawrap/refined_binning_{{sample}}/metawrap_" + str(config['metawrap_compl_thresh'])+ "_" + str(config['metawrap_contam_thresh']) + "_bins/bin.{i}.fa", 
#                       sample=config['sample'],
#                       i=range(1,6)), 
#    output:
#        directory("gtdbtk/{sample}/")
#    params:
#        extension     = "fa",
#        db_path     = config["gtdb_db"],
#        #mash_db     = addition to fast things up
#    threads: 16
#    resources:
#        mem_mb  = 200 * 1024, #200 GB
#        runtime = 4 * 60,
#    container: 
#        config["docker_gtdbtk"]
#    log: "logs/gtdbtk/{sample}.log"
#    shell:
#        """
#        export GTDBTK_DATA_PATH={params.db_path}
#        gtdbtk classify_wf --genome_dir {input.bins} --out_dir {output} --extension {params.extension} --cpus {threads} 2>&1 | tee {log}
#        """


#rule prokka:
#    input:
#        bins = expand("metawrap/refined_binning_{{sample}}/metawrap_" + str(config['metawrap_compl_thresh'])+ "_" + str(config['metawrap_contam_thresh']) + "_bins/bin.{i}.fa",
#                       sample=config['sample'],
#                       i=range(1,6)),
#    output:
#        directory("output_folder/prokka/{{sample}}/{bin}/")
#    params:
#        prefix = "{bin}",
#        mincontiglen=500,
#    threads: 8
#    resources:
#        mem_mb  = 200 * 1024, #200 GB
#        runtime = 4 * 60,
#    container:
#        config["docker_prokka"]
#    log: "logs/prokka/{sample}.log"
#    shell:
#        """
#        prokka --outdir {output} --prefix {params.prefix} --mincontiglen {params.mincontiglen} --cpus {threads} --force {input.bins} 2>&1 | tee {log}
#        """
 
rule gtdbtk_classify_wf:
    """
    Taxonomic classification of MAGs using GTDB-Tk.
    Input: refined MAG bins directory from MetaWRAP.
    Output: taxonomic classification results per sample.
    """
    input:
        #triggering 
        stats=f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats',
    output:
        directory("{sample}/gtdbtk/")
    params:
        bins_dir  = lambda wc: mag_dir(wc.sample),
        extension = "fa",
        db_path   = config["gtdb_db"],
        #mash_db   = config["gtdb_mash_db"],
    threads: 16
    resources:
        mem_mb  = 200 * 1024,
        runtime = 4 * 60,
    container:
        config["docker_gtdbtk"]
    log:
        "logs/gtdbtk/{sample}.log"
    shell:
        """
        export GTDBTK_DATA_PATH={params.db_path}
        gtdbtk classify_wf \
            --genome_dir {params.bins_dir} \
            --out_dir    {output}          \
            --extension  {params.extension} \
            --skip_ani_screen \
            --cpus       {threads}         \
            2>&1 | tee {log}
        """


rule prokka:
    """
    Functional annotation of a single MAG using Prokka.
    {sample} and {bin} wildcards — Snakemake runs this once per bin per sample.
    """
    input:
        bin = lambda wc: f"{mag_dir(wc.sample)}/{wc.bin}.fa",
    output:
        directory("{sample}/prokka/{bin}/")
    params:
        prefix       = "{bin}",
        mincontiglen = 500,
    threads: 8
    resources:
        mem_mb  = 16 * 1024,  # prokka is lightweight
        runtime = 2 * 60,
    container:
        config["docker_prokka"]
    log:
        "logs/prokka/{sample}/{bin}.log"
    shell:
        """
        prokka --outdir {output} --prefix {params.prefix} --metagenome --mincontiglen {params.mincontiglen} --cpus {threads} --force {input.bin} 2>&1 | tee {log}
        """


rule prokka_all_bins:
    """Aggregate all per-bin Prokka runs for a sample into a single done flag."""
    input:
        get_prokka_outputs
    output:
        touch("{sample}/prokka.done")

rule cleanup_rawbinning:
    input:
        # only runs after refinement is complete
        refined=f'{{sample}}/refined_binning/metawrap_{config["metawrap_compl_thresh"]}_{config["metawrap_contam_thresh"]}_bins.stats',
    output:
        touch("{sample}/logs/cleanup_rawbinning.done")
    params:
        rawdir="metawrap/rawbinning_{sample}/"
    shell:
        """
        rm -rf {params.rawdir}
        """
