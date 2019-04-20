#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType

from Bio import SeqIO


def parse_cds(file, key_id="protein_id", key_desc="product"):
	for record in SeqIO.parse(file, "genbank"):
		for feature in record.features:
			if feature.type == "CDS":
				cds = feature.extract(record)
				cds.id = feature.qualifiers[key_id][0]
				cds.description = feature.qualifiers[key_desc][0]
				yield cds


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
	parser.add_argument(
		"-key-id", "--key-id", default="protein_id"
	)
	parser.add_argument(
		"-key-description", "--key-description", dest="key_desc", default="product"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])
	with args.file as file:
		SeqIO.write(parse_cds(file, args.key_id, args.key_desc), sys.stdout, "fasta")


if __name__ == "__main__":
	sys.exit(main(sys.argv))
