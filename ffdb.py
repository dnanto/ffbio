#!/usr/bin/env python3

import sqlite3
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime
from pathlib import Path
from signal import signal, SIGPIPE, SIG_DFL
from subprocess import PIPE, Popen

from Bio import Entrez, SeqIO
from Bio.bgzf import BgzfWriter

formats = {
	"fas": "fasta",
	"faa": "fasta",
	"fna": "fasta",
	"fasta": "fasta",
	"gb": "gb",
	"gbk": "gb",
	"genbank": "gb"
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


def filepaths(path_idx):
	filenames = []

	if path_idx.exists():
		with sqlite3.connect(path_idx) as conn:
			rows = conn.execute("SELECT name FROM file_data ORDER BY file_number").fetchall()
			filenames = [path_idx.parent / row[0] for row in rows]

	return filenames


def reindex(path_idx, filenames, fmt, db, term, mdat, ext):
	path_idx_tmp = path_idx.with_suffix(".tmp")
	path_idx_tmp.exists() and path_idx_tmp.unlink()
	SeqIO.index_db(str(path_idx_tmp), filenames=filenames, format=fmt)
	with sqlite3.connect(path_idx_tmp) as conn:
		conn.execute(
			"INSERT INTO meta_data VALUES ('db', ?), ('term', ?), ('mdat', ?), ('ext', ?)",
			(db, term, mdat, ext)
		)
	path_idx_tmp.replace(path_idx)


def parse_argv(argv):
	parser = ArgumentParser(
		description="update a repository of indexed sequence files",
		formatter_class=ArgumentDefaultsHelpFormatter
	)

	parser.add_argument(
		"repo",
		type=Path,
		help="the target file to create/update"
	)
	parser.add_argument(
		"-fmt", "--fmt", "-format", "--format",
		default="fasta",
		help="the sequence file format"
	)
	parser.add_argument(
		"-db", "--db", "-database", "--database",
		default="nuccore",
		help="the NCBI database"
	)
	parser.add_argument(
		"-term", "--term",
		help="the NCBI query term"
	)
	parser.add_argument(
		"-no-mdat", "--no-mdat",
		action="store_true",
		help="flag to disable adding last modified date to the query"
	)
	parser.add_argument(
		"-email", "--email",
		default="",
		help="the e-mail to identify yourself to NCBI (for politeness reasons)"
	)
	parser.add_argument(
		"-post-size", "--post-size",
		type=int,
		default=1000,
		help="the number of records to post at a time"
	)
	parser.add_argument(
		"-cache", "--cache",
		nargs="+",
		help="the cache of sequence files to search prior to querying NCBI"
	)
	parser.add_argument(
		"-cache-fmt", "--cache-fmt",
		default="fasta",
		help="the sequence file format of the cache"
	)
	parser.add_argument(
		"-redo", "--redo",
		action="store_true",
		help="the flag to delete everything and redo"

	)

	args = parser.parse_args(argv)

	return args


def main(argv):
	args = parse_argv(argv[1:])

	Entrez.email = args.email

	repo = args.repo
	ext = args.fmt
	fmt = formats[ext]

	db = args.db
	term = args.term

	path_idx = repo.with_suffix(".idx")

	# load cache
	cache = {}
	if args.cache:
		cache = SeqIO.index_db(":memory:", filenames=args.cache, format=formats[args.cache_fmt])

	# redo or load metadata
	meta, files, accs = {}, [], set()
	if args.redo:
		for path in filepaths(path_idx):
			print("unlink: ", str(path))
			path.unlink()
		if path_idx.exists():
			print("unlink: ", str(path_idx))
			path_idx.unlink()
	elif path_idx.exists():
		with sqlite3.connect(path_idx) as conn:
			meta = dict(conn.execute("SELECT * FROM meta_data").fetchall())
			accs = {row[0] for row in conn.execute("SELECT key FROM offset_data")}

	db = meta.get("db", db)
	fmt = meta.get("format", fmt)
	ext = meta.get("ext", ext)
	mdat = None if args.no_mdat else meta.get("mdat")
	term = meta.get("term", term)
	term = f"{term} AND {mdat}:3000[MDAT]" if mdat else term

	print("term:", term, file=sys.stderr)

	# get the number of remote accessions
	mdat = datetime.now().strftime("%Y/%m/%d")
	with Entrez.esearch(db=db, term=term, idtype="acc") as file:
		record = Entrez.read(file)

	count = int(record["Count"])
	print("count: ", count, file=sys.stderr)

	# get the remote accessions
	with Entrez.esearch(db=db, term=term, idtype="acc", retmax=count) as file:
		record = Entrez.read(file)
		accs = set(record["IdList"]) - accs

	# accession diff length
	print("new: ", len(accs), file=sys.stderr)

	keys = accs & set(cache)
	accs -= keys
	print("cached:", len(keys), file=sys.stderr)
	print("download:", len(accs), file=sys.stderr)

	# update from cache
	if keys:
		filenames = list(map(str, filepaths(path_idx)))
		filenames += [str(path_idx.with_suffix(f".{len(filenames)}.{ext}.bgz"))]

		with BgzfWriter(filenames[-1]) as file:
			SeqIO.write((cache[key] for key in cache), file, fmt)

		reindex(path_idx, filenames, fmt, db, term, mdat, ext)

	if accs:
		# cat file | epost -db db | efetch -format fmt > target
		for batch in batchify(accs, args.post_size):
			print("download:", *batch[:5], "...", file=sys.stderr)
			kwargs = dict(stdin=PIPE, stdout=PIPE, universal_newlines=True)
			cmd1, cmd2 = ["epost", "-db", db], ["efetch", "-format", fmt]

			filenames = list(map(str, filepaths(path_idx)))
			filenames += [str(path_idx.with_suffix(f".{len(filenames)}.{ext}.bgz"))]

			with BgzfWriter(filenames[-1]) as file:
				with Popen(cmd1, **kwargs) as pipe1, Popen(cmd2, **kwargs) as pipe2:
					stdout, stderr = pipe1.communicate("\n".join(batch))
					stdout, stderr = pipe2.communicate(stdout)
					print(stdout, file=file)

			reindex(path_idx, filenames, fmt, db, term, mdat, ext)

	return 0


if __name__ == "__main__":
	signal(SIGPIPE, SIG_DFL)
	sys.exit(main(sys.argv))
