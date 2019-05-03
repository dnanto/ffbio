#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import OrderedDict
from signal import signal, SIGPIPE, SIG_DFL

import numpy as np
from Bio import SeqIO


def parse_argv(argv):
	parser = ArgumentParser(
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file", type=FileType()
	)
	parser.add_argument(
		"-fmt", "--fmt", "-format", "--format", default="fasta"
	)
	parser.add_argument(
		"-separator", "--separator", dest="sep", default="\t"
	)
	parser.add_argument(
		"-summary", "--summary", "-stats", "--stats", action="store_true"
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
