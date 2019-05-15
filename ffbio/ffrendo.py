#!/usr/bin/env python3


import sys
from argparse import ArgumentDefaultsHelpFormatter
from argparse import ArgumentParser, FileType
from collections import OrderedDict
from itertools import combinations
from signal import signal, SIGPIPE, SIG_DFL

from Bio import Restriction

from ffbio.ffparse import ffsniff, ffparse


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
	parser.add_argument(
		"-fragments", "--fragments",
		dest="frag",
		action="store_true",
		help="the flag to output fragments"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	# load enzymes
	enzymes = [getattr(Restriction, key) for key in args.enzymes]

	data = OrderedDict()
	with args.file as file:
		sample, fmt = ffsniff(file)
		for rec in ffparse(file, fmt, sample=sample):
			is_linear = rec.annotations.get("topology", not args.circular)
			data[rec.id] = {len(seq) for enz in enzymes for seq in enz.catalyze(rec.seq, is_linear)}

	for key1, key2 in combinations(data, 2):
		n = len(data[key1] & data[key2])
		d = 1 - (n / (len(data[key1]) + len(data[key2]) - n))
		print(key1, key2, d, sep="\t")

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
