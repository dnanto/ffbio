#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import defaultdict, OrderedDict

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid


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
		"-echo", "--echo",
		action="store_true",
		help="the flag to echo the sequence data to stdout"
	)
	parser.add_argument(
		"-prefix", "--prefix",
		help="the path prefix for the output files"
	)
	parser.add_argument(
		"-fmt-in", "--fmt-in", default="fasta"
	)
	parser.add_argument(
		"-fmt-out", "--fmt-out", default="fasta"
	)

	args = parser.parse_args(argv)
	args.prefix = args.prefix or args.file.name

	return args


def main(argv):
	args = parse_argv(argv[1:])

	# aggregate sequences by unique hash
	rec = defaultdict(list)
	with args.file as file:
		for record in SeqIO.parse(file, args.fmt_in):
			key = seguid(record.seq)
			rec[key].append(record)

	# map unique sequence index to unique hash
	idx_to_key = OrderedDict()
	with open(args.prefix + ".tsv", "w") as file:
		print("idx", "key", "id", "description", sep="\t", file=file)
		idx_width = len(str(len(rec)))
		for idx, ele in enumerate(rec.items(), start=1):
			key, val = ele
			idx_to_key[idx] = key
			for record in val:
				idx_val = "%0*d" % (idx_width, idx)
				print(idx_val, key, record.id, record.description, sep="\t", file=file)
				record.id, record.description = idx_val, ""

	# output unique sequences
	file = sys.stdout if args.echo else args.prefix + ".fna"
	SeqIO.write((rec[val][0] for val in idx_to_key.values()), file, args.fmt_out)

	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))
