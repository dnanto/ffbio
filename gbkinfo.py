#!/usr/bin/env python3

import sys
from Bio import SeqIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from signal import signal, SIGPIPE, SIG_DFL


def itemgetter(obj, *args):
	yield from (j.join(obj.get(k, d)) for k, d, j in args)


def parse_keys(keys):
	for key in keys:
		tokens = key.split(":", maxsplit=3)
		yield tokens + [""] * (3 - len(tokens))


def parse_argv(argv):
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		"file", type=FileType()
	)
	parser.add_argument(
		"keys", nargs="+"
	)
	parser.add_argument(
		"-fields", "-fields", "-header", "--header", dest="fields", action="store_true"
	)
	parser.add_argument(
		"-separator", "--separator", dest="sep", default="\t"
	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	keys = list(parse_keys(args.keys))

	with args.file as file:
		args.fields and print("id", *(key[0] for key in keys), sep=args.sep)
		for record in SeqIO.parse(file, "genbank"):
			feats = getattr(record, "features", [{}])
			quals = getattr(feats[0], "qualifiers", {})
			print(record.id, *itemgetter(quals, *keys), sep=args.sep)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
