#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from subprocess import PIPE, Popen

from Bio import Entrez
from Bio import SeqIO

ext_to_fmt = {
	"fna": "fasta",
	"gbk": "gb"
}


def batchify(entries, size=10):
	batch = []
	for i, e in enumerate(entries, start=1):
		batch.append(e)
		if i % size == 0:
			yield batch
			batch = []

	if batch:
		yield batch


def parse_argv(argv):
	parser = ArgumentParser(
		description="database auto util",
		formatter_class=ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		"target",
		type=Path,
		help="the target repor to create/update"
	)
	parser.add_argument(
		"term",
		help="the term"
	)
	parser.add_argument(
		"-db", "--db", "-database", "--database",
		default="nuccore",
		help="the database"
	)
	parser.add_argument(
		"-format", "--format",
		default="genbank",
		help="the sequence file format"
	)
	parser.add_argument(
		"-post_size", "--post-size",
		type=int,
		default=1000,
		help="the number of records to download at a time"
	)
	parser.add_argument(
		"-fetch-size", "--fetch-size",
		type=int,
		default=100,
		help="the number of records to download at a time"
	)
	parser.add_argument(
		"-email", "--email",
		default="",
		help="the e-mail"
	)
	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	Entrez.email = args.email

	target = args.target
	ext = target.name.split(".")[-1]
	fmt = ext_to_fmt.get(ext, ext)
	target.exists() or target.touch(exist_ok=True)
	accs = {rec.id for rec in SeqIO.parse(target, fmt)}

	with Entrez.esearch(db=args.db, term=args.term, idtype="acc") as file:
		record = Entrez.read(file)

	count = int(record["Count"])
	print("count: ", count, file=sys.stderr)
	with Entrez.esearch(db=args.db, term=args.term, idtype="acc", retmax=count) as file:
		record = Entrez.read(file)
		accs = set(record["IdList"]) - accs

	print("new: ", len(accs), file=sys.stderr)

	if accs:
		for batch in batchify(accs, args.post_size):
			print("download:", *batch[:5], "...", file=sys.stderr)
			kwargs = dict(stdin=PIPE, stdout=PIPE, universal_newlines=True)
			cmd1, cmd2 = ["epost", "-db", args.db], ["efetch", "-format", fmt]
			with target.open("a") as file:
				with Popen(cmd1, **kwargs) as pipe1:
					stdout, stderr = pipe1.communicate("\n".join(batch))
					with Popen(cmd2, **kwargs) as pipe2:
						stdout, stderr = pipe2.communicate(stdout)
						print(stdout, file=file)

	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))
