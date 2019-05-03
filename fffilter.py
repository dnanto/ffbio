#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def parse_argv(argv):
	parser = ArgumentParser(
		description="filter sequence records by header and/or length",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file",
		type=FileType(),
		help="the sequence file"
	)
	parser.add_argument(
		"-fmt", "--fmt", "-format", "--format",
		default="fasta",
		help="the sequence file format (input)"
	)
	parser.add_argument(
		"-fmt-o", "--fmt-o",
		default="fasta",
		help="the sequence file format (output)"
	)
	parser.add_argument(
		"-pattern",
		help="the regex pattern to search headers"
	)
	parser.add_argument(
		"-length",
		help="the sequence length"
	)
	parser.add_argument(
		"-percentage", "--percentage",
		action="store_true",
		help="the flag to use percentage bounds"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	rng = args.length
	if rng:
		rng = rng.split(":", maxsplit=3)
		rng += [0] * (3 - len(rng))
		base = int(rng[0])
		if args.percentage:
			rng = range(
				int(base - base * float(rng[2]) / 100),
				int(base + base * float(rng[1]) / 100)
			)
		else:
			rng = range(base - int(rng[1]), base + int(rng[2]))

	with args.file as file:
		records = SeqIO.parse(file, args.fmt)
		if rng:
			records = (rec for rec in records if len(rec) in rng)
		if args.pattern:
			records = (rec for rec in records if re.search(args.pattern, rec.description))
		SeqIO.write(records, sys.stdout, args.fmt_o)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
