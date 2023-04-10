#! /usr/bin/env python3

import argparse
import csv

# parses command-line arguments with DESeq csv file, original gff file, counts matrix, and output file
def parse_args():
    parser = argparse.ArgumentParser(description="Filters gff file based on filtered fasta file")
    parser.add_argument("deseq", help="csv file in standard deseq format")
    parser.add_argument("--gff", help="gff file containing functional annotations")
    parser.add_argument("--counts", help="counts matrix containing read counts")
    parser.add_argument("--out", "-o", default="deseq_summary", help="prefix for output files")
    return parser.parse_args()

# parse deseq csv
def parse_deseq_csv(deseq_csv):
    # [{"": "ATP6", "baseMean": 14, "log2FoldChange": -1.44, ...}, ...]
    matrix_data = []

    with open(deseq_csv) as matrix:
        reader = csv.DictReader(matrix)
        for line in reader:
            matrix_data.append(line)

    return matrix_data

# parse the count matrix
def parse_count_matrix(count_matrix):
    # {"ATP6": {"": "ATP6", "EV_1": 14, ...}, ...}
    matrix_data = {}

    with open(count_matrix) as matrix:
        reader = csv.DictReader(matrix, dialect="excel-tab")
        for line in reader:
            matrix_data.setdefault(line[""], line)

    return matrix_data

# parse the gff file to get functional annotations of all genes
def parse_gff(gff_file):
    gff_dict = {}

    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#"):
                continue

            attributes = line.rstrip().split("\t")[-1]

            attribute_list = (";").join(attributes.split("=")).split(";")

            # look into attributes to find locus and product
            try:
                locus = attribute_list.index("locus_tag")+1
            except ValueError:
                continue

            try:
                product = attribute_list.index("product")+1
            except ValueError:
                if "pseudogene" in attribute_list:
                    product = attribute_list.index("pseudogene")
                else:
                    continue
            
            gff_dict.setdefault(attribute_list[locus], attribute_list[product])

    return gff_dict

# write out new csv with data
def write_summary(deseq_data, count_data, gff_data, out_prefix):
    # deseq_data: [{"": "ATP6", "baseMean": 14, "log2FoldChange": -1.44, ...}, ...]
    # count_data: {"ATP6": {"": "ATP6", "EV_1": 14, ...}, ...}
    # gff_data: {"ATP6": "large subunit of ATP synthase"}

    with open("{0}.csv".format(out_prefix), "w") as csvfile:
        fieldnames = []
        for feat_dict in deseq_data:
            output_dict = {}

            # extract feature and annotation
            feature = feat_dict.pop("")
            output_dict["feature"] = feature
            try:
                output_dict["annotation"] = gff_data[feature]
            except KeyError:
                output_dict["annotation"] = ""

            # extract counts per sample
            count_dict = count_data[feature]
            count_dict.pop("")
            output_dict = dict(output_dict, **count_dict, **feat_dict)

            # write header first time only
            if len(fieldnames) == 0:
                fieldnames = list(output_dict.keys())
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()

            writer.writerow(output_dict)

def main():
    args = parse_args()
    deseq_data = parse_deseq_csv(args.deseq)
    count_data = parse_count_matrix(args.counts)
    gff_data = parse_gff(args.gff)
    write_summary(deseq_data, count_data, gff_data, args.out)

main()