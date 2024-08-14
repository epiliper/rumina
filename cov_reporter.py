import pandas as pd
import os
import warnings
import pysam

warnings.simplefilter(action="ignore", category=FutureWarning)

COLUMNS = [
    "num_reads_input_file",
    "num_raw_umis",
    "num_total_groups",
    "num_passing_groups",
    "least_reads_group",
    "least_reads_per_group",
    "most_reads_group",
    "most_reads_per_group",
]


def generate_report(original_file, final_file):
    work_path = os.path.dirname(final_file)
    minmax_file = os.path.join(work_path, "minmax.txt")

    df = pd.read_csv(
        minmax_file,
        sep="\t",
        names=[
            "min_groups",
            "mins",
            "max_groups",
            "maxes",
            "num_passing_groups",
            "num_total_groups",
            "num_umis",
        ],
    )

    report = pd.DataFrame(columns=COLUMNS)
    report["least_reads_per_group"] = [int(df["mins"].min())]

    report["least_reads_group"] = [df.iloc[df["mins"].idxmin()].min_groups]

    report["most_reads_per_group"] = [int(df["maxes"].max())]

    report["most_reads_group"] = [df.iloc[df["maxes"].idxmax()].max_groups]

    report["num_passing_groups"] = [df["num_passing_groups"].sum()]

    report["num_total_groups"] = [df["num_total_groups"].sum()]

    report["num_raw_umis"] = [df["num_umis"].sum()]

    report["num_reads_input_file"] = [
        pysam.AlignmentFile(original_file).count(until_eof=True)
    ]

    report["num_reads_output_file"] = [
        pysam.AlignmentFile(final_file).count(until_eof=True)
    ]

    print(f"written {report["num_reads_output_file"].values[0]} reads...")

    try:
        cov_stats = pysam.coverage(final_file).split("\n")
        cov_stats = zip(cov_stats[0].split("\t"), cov_stats[1].split("\t"))

        for col_val in cov_stats:
            report[col_val[0]] = [col_val[1]]

    except pysam.utils.SamtoolsError:
        print(
            "coverage and depth reporting failed. Output BAM could not be read by pysam."
        )

    save_dir = os.path.dirname(final_file)
    report["query_name"] = os.path.basename(original_file).split(".bam")[0]

    csv_name = os.path.join(save_dir, final_file.split(".bam")[0] + "_coverage.tsv")
    report.to_csv(csv_name, sep="\t", index=None)

    # IMPORTANT:
    # for accurate reporting of split files,
    # rumina will append to minmax file,
    # delete it after each sample
    os.remove(minmax_file)


def summarize_coverage(work_dir):
    if not os.path.isabs(work_dir):
        work_dir = os.path.abspath(work_dir)

    total_df = pd.DataFrame(columns=COLUMNS)

    final_file = os.path.join(work_dir, "COVERAGE_REPORT.csv")

    print(f"Combining coverage reports...\nCoverage csv name: {final_file}")

    for file in os.listdir(work_dir):
        if file.endswith("_coverage.tsv"):
            read_file = pd.read_csv(
                os.path.join(work_dir, file), sep="\t", header=0, index_col=False
            )
    total_df = pd.concat([total_df, read_file])

    total_df.to_csv(final_file)

    for file in os.listdir(work_dir):
        if file.endswith("_coverage.tsv"):
            os.remove(os.path.join(work_dir, file))
