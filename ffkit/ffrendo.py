#!/usr/bin/env python3


import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import Restriction

from ffkit.ffparse import ffsniff, ffparse


def parse_argv(argv):
	parser = ArgumentParser(
		description="compute restriction fragment length polymorphism profile",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file",
		type=FileType(),
		help="the sequence file"
	)
	parser.add_argument(
		"enzymes",
		nargs="+",
		help="the restriction enzyme names"
	)
	parser.add_argument(
		"-circular", "--circular",
		action="store_true",
		help="the flag to set default topology to circular if not specified in annotation"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	# load enzymes
	enzymes = [getattr(Restriction, key) for key in args.enzymes]

	# aggregate sequences by unique hash
	with args.file as file:
		sample, fmt = ffsniff(file)
		for record in ffparse(file, fmt, sample=sample):
			for enzyme in enzymes:
				is_linear = record.annotations.get("topology", not args.circular)
				for seq in enzyme.catalyze(record.seq, is_linear):
					print(record.id, str(enzyme), len(seq), sep="\t")

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
