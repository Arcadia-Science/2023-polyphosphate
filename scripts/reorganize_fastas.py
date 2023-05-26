#! /usr/bin/python3

from Bio import SeqIO
import os
import argparse


def main():
    # command line arguments
    parser = argparse.ArgumentParser(description='Rename FASTA headers for files in a directory and concatenate to single FASTAs')
    parser.add_argument('input_dir', help='Directory containing input FASTA files')
    parser.add_argument('output_file', help='Output concatenated FASTA file')
    args = parser.parse_args()

    # dictionary to hold records
    new_records = {}

    # loop through files in a directory
    for filename in os.listdir(args.input_dir):
        # skip non-fasta files
        if not filename.endswith(".fasta"):
            continue
        # remove directory path and extension so just have accession name
        accession_name = os.path.splitext(os.path.basename(filename))[0]
        # open fasta file and modify the headers
        with open(os.path.join(args.input_dir, filename), "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                record.id = f"{accession_name}"
                record.description=""
                new_records[record.id] = record

    # write to the new concatenated FASTA file
    with open(args.output_file, "w") as outfile:
        SeqIO.write(new_records.values(), outfile, "fasta")


if __name__ == '__main__':
    main()
