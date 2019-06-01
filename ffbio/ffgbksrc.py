#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO


def itemgetter(obj, *args):
	yield from (j.join(obj.get(k, d)) for k, d, j in args)


def parse_keys(keys):
	for key in keys:
		tokens = key.split(":", maxsplit=3)
		yield tokens + [""] * (3 - len(tokens))


def parse_argv(argv):
	parser = ArgumentParser(
		description="output a table of source feature metadata from a GenBank file",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"file",
		type=FileType(),
		help="the sequence file"
	)
	parser.add_argument(
		"keys",
		nargs="+",
		help="the qualifier keys"
	)
	parser.add_argument(
		"-fields", "-fields", "-header", "--header",
		dest="fields",
		action="store_true",
		help="the flag to output a header"
	)
	parser.add_argument(
		"-separator", "--separator",
		dest="sep",
		default="\t",
		help="the table delimiter, the default is a tab character"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	keys = list(parse_keys(args.keys))

	with args.file as file:
		args.fields and print("id", "description", *(key[0] for key in keys), sep=args.sep)
		for record in SeqIO.parse(file, "genbank"):
			feats = getattr(record, "features", [{}])
			quals = getattr(feats[0], "qualifiers", {})
			print(record.id, record.description, *itemgetter(quals, *keys), sep=args.sep)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
