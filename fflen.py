#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import OrderedDict
from signal import signal, SIGPIPE, SIG_DFL

import numpy as np
from Bio import SeqIO


def parse_argv(argv):
	parser = ArgumentParser(
		description="compute the length of each record in the file",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file", type=FileType(),
		help="the sequence file"
	)
	parser.add_argument(
		"-fmt", "--fmt", "-format", "--format",
		default="fasta",
		help="the sequence file format"
	)
	parser.add_argument(
		"-separator", "--separator",
		dest="sep",
		default="\t",
		help="the table delimiter, default is the tab character"
	)
	parser.add_argument(
		"-summary", "--summary", "-stats", "--stats",
		action="store_true",
		help="the flag to compute basic summary statistics"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	with args.file as file:
		records = SeqIO.parse(file, args.fmt)
		if args.summary:
			lens = np.array([len(rec) for rec in records])
			results = OrderedDict(
				(
					("mean", str(lens.mean())),
					("median", str(np.median(lens))),
					("std", str(lens.std())),
					("min", str(lens.min())),
					("max", str(lens.max()))
				)
			)
			print(*(args.sep.join(item) for item in results.items()), sep="\n")
		else:
			for record in records:
				print(record.id, record.description, len(record), sep=args.sep)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
