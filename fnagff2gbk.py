#!/usr/bin/env python3

import fileinput
import os.path
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation


def parse_argv(argv):
	parser = ArgumentParser(
		description="compute the unique set of sequences",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"-fna",
		nargs="+",
		required=True,
		help="the sequence files"
	)
	parser.add_argument(
		"-gff",
		nargs="+",
		required=True,
		help="the gff file"
	)
	parser.add_argument(
		"-out", "--out",
		type=Path,
		default=Path()
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	records = SeqIO.index_db(
		":memory:",
		filenames=args.fna,
		format="fasta",
		alphabet=generic_dna
	)

	os.makedirs(args.out, exist_ok=True)

	features = defaultdict(list)
	with fileinput.input(args.gff) as file:
		for line in file:
			tokens = line.rstrip().split("\t")
			seqid, source, type, start, end, score, strand, phase, atts = tokens
			atts = dict(ele.split("=") for ele in atts.split(";"))
			atts = {key: val.strip() for key, val in atts.items()}
			location = FeatureLocation(start=int(start), end=int(end), strand=[-1, 1][strand == "+"])
			feature = SeqFeature(location, type=type, qualifiers=atts)
			features[seqid].append(feature)

	for key, val in records.items():
		path = args.out.joinpath(key).with_suffix(".gbk")
		with path.open("w") as file:
			val.name = key
			val.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
			val.features = sorted(features[key], key=lambda x: x.location.start)
			print(path)
			SeqIO.write(val, file, "genbank")


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
