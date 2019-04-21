#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def parse_cds(file):
	for record in SeqIO.parse(file, "genbank"):
		for feature in record.features:
			if feature.type == "CDS":
				try:
					protein_id = feature.qualifiers["protein_id"][0]
					product = feature.qualifiers.get("product", ["n/a"])[0]
					cds = feature.extract(record)
					cds.id = f"lcl|{protein_id}"
					cds.description = f"{record.id}|{product}|{feature.location}"
					yield cds
				except:
					pass


def parse_argv(argv):
	parser = ArgumentParser(
		description="compute the unique set of sequences",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file",
		type=FileType(),
		help="the sequence file"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])
	with args.file as file:
		SeqIO.write(parse_cds(file), sys.stdout, "fasta")


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
