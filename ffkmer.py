#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def kmerize(seq, k=2, step=1):
	yield from (seq[i:i + k] for i in range(0, len(seq) - k + 1, step))


def parse_argv(argv):
	parser = ArgumentParser(
		description="calculate the set of k-mers",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file",
		type=FileType(),
		help="the sequence file"
	)
	parser.add_argument(
		"-fmt", "--fmt", "-format", "--format",
		dest="fmt",
		default="fasta",
		help="the sequence file format"
	)
	parser.add_argument(
		"size",
		type=int,
		help="the k-mer size"
	)
	parser.add_argument(
		"-step",
		type=int,
		default=1,
		help="the next position to scan in the sequence"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	with args.file as file:
		for record in SeqIO.parse(file, args.fmt):
			for idx, ele in enumerate(kmerize(str(record.seq), args.size, args.step)):
				print(f">{record.id}-{idx}", ele, sep="\n")

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
