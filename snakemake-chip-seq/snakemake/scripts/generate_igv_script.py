import os
import random
import math

import pandas as pd
import pyBigWig

from typing import Optional, List

random.seed(123)

batch_script = []
batch_script.append(f"genome {snakemake.params.genome_label}")
if "bam" in snakemake.input.keys():
    batch_script.append("preference SAM.SHOW_ALIGNMENT_TRACK false")

sample_table = pd.DataFrame(snakemake.params.sample_table)

filenames = sample_table["filename"].replace(
    r"\/?([^\/]+)\.f(ast)?q(?:\.gz)?$", r"\1", regex=True
)
genotypes = set(sample_table["genotype"].tolist())

genotype_colors = {
    genotype: (
        random.randint(1, 256),
        random.randint(1, 256),
        random.randint(1, 256),
    )
    for genotype in genotypes
}

external_tracks_color = (
    random.randint(1, 256),
    random.randint(1, 256),
    random.randint(1, 256),
)

with open(snakemake.input.genome) as fin:
    chromsizes = {}
    for line in fin:
        line = line.strip().split("\t")
        chromsizes[line[0]] = int(line[1])

coverage = {}
if "coverage" in snakemake.input.keys() and "bam" in snakemake.input.keys():
    with open(snakemake.input.coverage) as fin:
        for line in fin:
            try:
                line = line.strip().split("\t")
                key = f"{line[0]}:{line[1]}-{line[2]}"
                cov = [int(x) for x in line[4:]]
                cov = max(cov) + 2
                coverage[key] = cov
            except IndexError as err:
                print(line)
                raise err


def get_sample(p, ext):
    filename = os.path.basename(p).replace(ext, "")
    df = sample_table[filenames == filename]
    return df


def get_track_name(p, ext):
    df = get_sample(p, ext)
    track_name = df["sample_name"].iloc[0]
    return track_name


def get_external_track_name(p, ext):
    filename = os.path.basename(p).replace(ext, "")
    return filename


def get_track_color(p, ext):
    df = get_sample(p, ext)
    genotype = df["genotype"].iloc[0]
    track_color = ",".join([str(x) for x in genotype_colors[genotype]])
    return track_color


def calc_sloped_region(chrom, start, end):
    start = start - 1000 if start - 1000 > 1 else start
    end = end + 1000 if end + 1000 < chromsizes[chrom] else end
    return chrom, start, end


def get_formatted_region(chrom, start, end):
    chrom, start, end = calc_sloped_region(chrom, start, end)
    formatted_region = f"{chrom}:{start}-{end}"
    return formatted_region


def get_interval_range(chrom, start, end, bw_handlers: List = [], extra=2):
    if "bam" in snakemake.input.keys():
        key = f"{chrom}:{start}-{end}"
        ret = coverage[key]
    else:
        chrom, start, end = calc_sloped_region(chrom, start, end)
        bw_max_values = []
        for bw in bw_handlers:
            max_value = bw.stats(chrom, start, end, type="max", exact=True)
            bw_max_values += max_value
        ret = max(bw_max_values)
        ret = math.ceil(ret) + extra
    return ret


track_names = []
bw_handlers = []
if "bam" in snakemake.input.keys():
    for bam in snakemake.input.bam:
        track_name = get_track_name(bam, ".bam")
        track_color = get_track_color(bam, ".bam")

        batch_script.append(f"load {bam} name={track_name}")
        batch_script.append(f"setColor {track_color} {track_name}")
else:
    for bw in snakemake.input.bigwig:
        track_name = get_track_name(bw, ".bw")
        track_color = get_track_color(bw, ".bw")

        batch_script.append(f"load {bw} name={track_name}")
        batch_script.append(f"setColor {track_color} {track_name}")
        batch_script.append(f"setLogScale false {track_name}")

        bw_handlers.append(pyBigWig.open(bw))
        track_names.append(track_name)

    if "external_dataset" in snakemake.input.keys():
        external_bw_handlers = []
        external_track_names = []
        for bw in snakemake.input.external_dataset:
            track_name = get_external_track_name(bw, ".bw")
            track_color = external_tracks_color

            batch_script.append(f"load {bw} name={track_name}")
            batch_script.append(f"setColor {track_color} {track_name}")
            batch_script.append(f"setLogScale false {track_name}")

            external_bw_handlers.append(pyBigWig.open(bw))
            external_track_names.append(track_name)


for bed in snakemake.input.bed:
    track_name = get_track_name(bed, "_peaks.bed")
    track_color = get_track_color(bed, "_peaks.bed")

    batch_script.append(f"load {bed} name={track_name} peaks")
    batch_script.append(f"setColor {track_color} {track_name}")

batch_script.append(f"load {snakemake.input.calls} name=Overlaps")


with open(snakemake.input.calls) as fin:
    regions = fin.readlines()
    for picture_format in ["svg", "png"]:
        batch_script.append(
            f"snapshotDirectory {snakemake.params.snapshot_directory}/{picture_format}"
        )
        for region in regions:

            region = region.rstrip().split("\t")
            chrom, start, end, region_name, _ = region
            start = int(start)
            end = int(end)

            formatted_region = get_formatted_region(chrom, start, end)
            batch_script.append(f"goto {formatted_region}")

            interval_range = get_interval_range(chrom, start, end, bw_handlers)

            for track_name in track_names:
                batch_script.append(f"setDataRange 0,{interval_range:.2f} {track_name}")

            if "external_dataset" in snakemake.input.keys():
                external_tracks_data_range = get_interval_range(
                    chrom, start, end, external_bw_handlers
                )
                for track_name in external_track_names:
                    batch_script.append(
                        f"setDataRange 0,{external_tracks_data_range:.2f} {track_name}"
                    )

            batch_script.append(f"snapshot {region_name}.{picture_format}")

batch_script.append("exit")

with open(snakemake.output.batch, "w") as fout:
    fout.write("\n".join(batch_script))
