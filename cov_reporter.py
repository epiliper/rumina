import pandas as pd
import os
import warnings
import pysam
from concurrent.futures import ThreadPoolExecutor

warnings.simplefilter(action="ignore", category=FutureWarning)

COLUMNS = [
    "num_reads_input_file",
    "num_reads_output_file",
    "num_raw_umis",
    "num_total_groups",
    "num_passing_groups",
    "least_reads_group",
    "least_reads_per_group",
    "most_reads_group",
    "most_reads_per_group",
]


def get_coverage(input_file, region):
    region_df = pd.DataFrame()

    cov_stats = pysam.coverage("-r", region, input_file)
    cov_stats = cov_stats.split("\n")
    cov_stats = zip(cov_stats[0].split("\t"), cov_stats[1].split("\t"))

    for col_val in cov_stats:
        region_df[col_val[0]] = [col_val[1]]

    return region_df


def generate_report(original_file, final_file):
    work_path = os.path.dirname(final_file)
    minmax_file = os.path.join(work_path, "minmax.txt")

    df = pd.read_csv(
        minmax_file,
        sep="\t",
        names=COLUMNS,
    )

    group_report = pd.DataFrame(columns=COLUMNS)
    cov_report = pd.DataFrame()

    group_report["num_reads_input_file"] = [df["num_reads_input_file"].sum()]
    group_report["num_reads_output_file"] = [df["num_reads_output_file"].sum()]
    group_report["least_reads_per_group"] = [int(df["least_reads_per_group"].min())]
    group_report["least_reads_group"] = [
        df.iloc[df["least_reads_per_group"].idxmin()].least_reads_group
    ]
    group_report["most_reads_per_group"] = [int(df["most_reads_per_group"].max())]
    group_report["most_reads_group"] = [
        df.iloc[df["most_reads_per_group"].idxmax()].most_reads_group
    ]
    group_report["num_passing_groups"] = [df["num_passing_groups"].sum()]
    group_report["num_total_groups"] = [df["num_total_groups"].sum()]
    group_report["num_raw_umis"] = [df["num_raw_umis"].sum()]

    print(f"written {group_report["num_reads_output_file"].values[0]} reads...")

    with pysam.AlignmentFile(final_file, "rb") as outbam:
        references = outbam.references

    # use a thread for each region
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(get_coverage, final_file, region) for region in references
        ]
        for future in futures:
            cov_report = pd.concat([cov_report, future.result()], ignore_index=True)

    save_dir = os.path.dirname(final_file)
    group_report["query_name"] = os.path.basename(original_file).split(".bam")[0]

    group_tsv = os.path.join(save_dir, final_file.split(".bam")[0] + "_grouping.tsv")
    cov_tsv = os.path.join(save_dir, final_file.split(".bam")[0] + "_coverage.tsv")

    group_report.to_csv(group_tsv, sep="\t", index=None)
    cov_report.to_csv(cov_tsv, sep="\t", index=None)

    # IMPORTANT:
    # for accurate reporting of split files,
    # rumina will append to minmax file,
    # delete it after each sample
    os.remove(minmax_file)
