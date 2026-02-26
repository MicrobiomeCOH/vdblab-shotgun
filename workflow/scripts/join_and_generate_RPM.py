import sys
import subprocess
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def concat_tsv_files(file_list):
    """Read a list of TSV files and concatenate, keeping one header."""
    dfs = []
    for path in file_list:
        try:
            df = pd.read_csv(path, sep="\t")
            if not df.empty:
                dfs.append(df)
        except pd.errors.EmptyDataError:
            # touched-but-empty files from sparse batches are expected
            pass
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()


def count_num_reads_compressed_file(file_name):
    """Return the number of reads in a gzipped FASTQ file."""
    command = f"echo $(($(zcat {file_name} | wc -l)/4))"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    return float(result.stdout.strip())


def build_rpm_table(overview_df, substrate_df, cgc_df, coverage_path, num_reads):
    """
    Merge all CAZI tables with coverage and compute RPM.

    coverage_path is the output of:
        bedtools genomecov -d -5 -ibam sample.bam
    which produces 3 columns: gene_name, position, depth  (no header).
    """
    # ------------------------------------------------------------------
    # Coverage: sum per-position depth to get total counts per gene
    # ------------------------------------------------------------------
    counts_df = pd.read_csv(
        coverage_path,
        sep="\t",
        header=None,
        names=["bam_name", "position", "depth"],
    )
    counts_df = (
        counts_df.groupby("bam_name")
        .agg(counts=("depth", "sum"), length=("position", "max"))
        .reset_index()
    )

    # ------------------------------------------------------------------
    # Substrate: split the composite key_name column
    # ------------------------------------------------------------------
    if not substrate_df.empty:
        sub_df = substrate_df.copy()
        # Handle both named and unnamed substrate files
        if "k_name_parts" not in sub_df.columns:
            sub_df.columns = [
                "k_name_parts",
                "PULID",
                "substrate",
                "substrate_bitscore",
                "signature pairs",
                "dbCAN-sub substrate",
                "dbCAN-sub substrate score",
            ]
        sub_df["k_name"] = sub_df["k_name_parts"].map(lambda x: x.split("|")[0])
        sub_df["cgc"]    = sub_df["k_name_parts"].map(lambda x: x.split("|")[1])
        sub_df = sub_df[["k_name", "cgc", "substrate", "substrate_bitscore"]]
    else:
        sub_df = pd.DataFrame(
            columns=["k_name", "cgc", "substrate", "substrate_bitscore"]
        )

    # ------------------------------------------------------------------
    # CGC: key file linking k_names <-> bam_names
    # ------------------------------------------------------------------
    if not cgc_df.empty:
        cgc_faa_df = cgc_df.copy()
        # raw cgc.out has no header; rename if needed
        if cgc_faa_df.columns[0] == 0 or not isinstance(cgc_faa_df.columns[0], str):
            cgc_faa_df.columns = [
                "_0", "caz_type", "_2", "_3",
                "cgc", "k_name", "_6", "_7",
                "bam_name", "_9", "_10", "_11",
            ]
        cgc_faa_df = cgc_faa_df[["caz_type", "cgc", "k_name", "bam_name"]]
    else:
        cgc_faa_df = pd.DataFrame(
            columns=["caz_type", "cgc", "k_name", "bam_name"]
        )

    # ------------------------------------------------------------------
    # Merge everything together
    # ------------------------------------------------------------------
    caz_df = overview_df.merge(
        counts_df, how="left", left_on="Gene ID", right_on="bam_name"
    ).drop("bam_name", axis="columns", errors="ignore")

    caz_df = caz_df.merge(
        cgc_faa_df, how="left", left_on="Gene ID", right_on="bam_name"
    )

    caz_df = caz_df.merge(sub_df, how="left", on=["k_name", "cgc"])

    caz_df["RPM"] = ((10 ** 6) * caz_df["counts"]) / num_reads

    return caz_df


def read_cgc_raw(file_list):
    """
    cgc.out has no real header and uses '+' lines as separators.
    Read all batch files, skip comment lines, concatenate.
    """
    dfs = []
    col_names = [
        "_0", "caz_type", "_2", "_3",
        "cgc", "k_name", "_6", "_7",
        "bam_name", "_9", "_10", "_11",
    ]
    for path in file_list:
        try:
            df = pd.read_csv(
                path,
                sep="\t",
                comment="+",
                header=None,
                names=col_names,
            )
            if not df.empty:
                dfs.append(df)
        except pd.errors.EmptyDataError:
            pass
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame(columns=col_names)


def read_substrate_raw(file_list):
    """
    substrate.out has a header row but it is not always present in empty files.
    """
    col_names = [
        "k_name_parts",
        "PULID",
        "substrate",
        "substrate_bitscore",
        "signature pairs",
        "dbCAN-sub substrate",
        "dbCAN-sub substrate score",
    ]
    dfs = []
    for path in file_list:
        try:
            df = pd.read_csv(path, sep="\t", skiprows=1, header=None, names=col_names)
            if not df.empty:
                dfs.append(df)
        except pd.errors.EmptyDataError:
            pass
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame(columns=col_names)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(
    overview_files,
    substrate_files,
    cgc_files,
    coverage_path,
    r1_path,
    out_overview,
    out_substrate,
    out_cgc,
    out_rpm,
):
    # --- merge per-batch flat files ---
    overview_df  = concat_tsv_files(overview_files)
    substrate_df = read_substrate_raw(substrate_files)
    cgc_df       = read_cgc_raw(cgc_files)

    # --- write merged flat outputs ---
    overview_df.to_csv(out_overview, sep="\t", index=False)
    substrate_df.to_csv(out_substrate, sep="\t", index=False)
    cgc_df.to_csv(out_cgc, sep="\t", index=False)

    # --- compute RPM ---
    if overview_df.empty:
        print("WARNING: no CAZyme hits found across all batches; RPM file will be empty")
        pd.DataFrame().to_csv(out_rpm, sep="\t")
        return

    num_reads = count_num_reads_compressed_file(r1_path)
    print(f"Total reads in {r1_path}: {num_reads:.0f}")

    rpm_df = build_rpm_table(
        overview_df, substrate_df, cgc_df, coverage_path, num_reads
    )
    rpm_df.to_csv(out_rpm, sep="\t", index=False)
    print(f"Wrote RPM table with {len(rpm_df)} rows to {out_rpm}")


if __name__ == "__main__":
    if "snakemake" in globals():
        # Snakemake passes lists for multi-input rules
        overview_files  = list(snakemake.input.overview)
        substrate_files = list(snakemake.input.substrate)
        cgc_files       = list(snakemake.input.cgc)
        coverage_path   = str(snakemake.input.coverage)
        r1_path         = str(snakemake.input.r1)
        out_overview    = str(snakemake.output.overview)
        out_substrate   = str(snakemake.output.substrate)
        out_cgc         = str(snakemake.output.cgc)
        out_rpm         = str(snakemake.output.rpm)
    else:
        # CLI fallback: overview [overview ...] substrate [substrate ...]
        # Pass files as positional args in groups separated by '--'
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument("--overview",  nargs="+", required=True)
        p.add_argument("--substrate", nargs="+", required=True)
        p.add_argument("--cgc",       nargs="+", required=True)
        p.add_argument("--coverage",  required=True)
        p.add_argument("--r1",        required=True)
        p.add_argument("--out-overview",  required=True, dest="out_overview")
        p.add_argument("--out-substrate", required=True, dest="out_substrate")
        p.add_argument("--out-cgc",       required=True, dest="out_cgc")
        p.add_argument("--out-rpm",       required=True, dest="out_rpm")
        args = p.parse_args()
        overview_files  = args.overview
        substrate_files = args.substrate
        cgc_files       = args.cgc
        coverage_path   = args.coverage
        r1_path         = args.r1
        out_overview    = args.out_overview
        out_substrate   = args.out_substrate
        out_cgc         = args.out_cgc
        out_rpm         = args.out_rpm

    main(
        overview_files,
        substrate_files,
        cgc_files,
        coverage_path,
        r1_path,
        out_overview,
        out_substrate,
        out_cgc,
        out_rpm,
    )
