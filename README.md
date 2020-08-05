---
title: "README"
author: "dnanto"
date: "August 05, 2020"
output: 
  html_document: 
    keep_md: yes
    toc: yes
    toc_depth: 2
---



# Abstract

**ffbio** is a collection of scripts to work with flat-file sequence databases and help me complete my thesis.

Note: this repository is experimental...

Note: this is an Rmarkdown document and uses a mix of R and bash code chunks.

* The Pipfile lists the Python requirements.
* For convenience, add the script root directory to PATH.
* Also, it is convenient to know the NCBI [search field descriptions](https://www.ncbi.nlm.nih.gov/books/NBK49540/).

# Scripts

## ffdb.py

The ```ffdb.py``` script searches NCBI and downloads the entire result set. The destination folder stores the [BGZip](http://www.htslib.org/doc/bgzip.html)-compressed flat-files. The [index_db](https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html#index_db) method generates the [SQLite](https://www.sqlite.org/index.html) database in the same directory. The script also modifies the meta_data table to preserve NCBI search parameters.


```bash
ffdb.py -h
```

```
usage: ffdb.py [-h] [-db DB] [-term TERM] [-rettype RETTYPE] [-retmax RETMAX]
               [-xmdat] [-email EMAIL]
               repo

update a repository of indexed sequence files

positional arguments:
  repo              the target file to create/update

optional arguments:
  -h, --help        show this help message and exit
  -db DB            the NCBI database (default: nuccore)
  -term TERM        the NCBI query term (default: None)
  -rettype RETTYPE  the sequence file format (default: fasta)
  -retmax RETMAX    the records to post at a time (default: 1000)
  -xmdat            flag to disable modified date limit (default: False)
  -email EMAIL      the e-mail to identify yourself to NCBI (default: )
```

### Ex-1

Initialize a repository of indexed GenBank files. Search NCBI for all *Mimivirus* genomic DNA sequences.


```bash
ffdb.py data/315393 \
  -term 'txid315393[PORGN] AND biomol_genomic[PROP] NOT gbdiv_pat[PROP] NOT gbdiv_syn[PROP] NOT WGS[KYWD]' \
  -rettype fasta \
  -retmax 100
```

```
2020-08-05T11:49:48 ['/Users/dnanto/GitHub/ffbio/ffbio/ffdb.py', 'data/315393', '-term', 'txid315393[PORGN] AND biomol_genomic[PROP] NOT gbdiv_pat[PROP] NOT gbdiv_syn[PROP] NOT WGS[KYWD]', '-rettype', 'fasta', '-retmax', '100']
2020-08-05T11:49:48 txid315393[PORGN] AND biomol_genomic[PROP] NOT gbdiv_pat[PROP] NOT gbdiv_syn[PROP] NOT WGS[KYWD]
2020-08-05T11:49:50 count = 355
2020-08-05T11:49:57 000 - 100 028.17%
2020-08-05T11:50:07 100 - 200 056.34%
2020-08-05T11:50:13 200 - 300 084.51%
2020-08-05T11:50:17 300 - 355 100.00%
2020-08-05T11:50:17 data/315393-1.fasta.bgz.tmp -> data/315393-1.fasta.bgz
2020-08-05T11:50:17 data/315393-2.fasta.bgz.tmp -> data/315393-2.fasta.bgz
2020-08-05T11:50:17 data/315393-3.fasta.bgz.tmp -> data/315393-3.fasta.bgz
2020-08-05T11:50:17 data/315393-4.fasta.bgz.tmp -> data/315393-4.fasta.bgz
index...
```

List the results in the directory.


```bash
du -sh data/315393.*
```

```
 40K	data/315393.db
4.0K	data/315393.log
```

### Ex-2

Update the same repository.


```bash
ffdb.py data/315393
```

```
2020-08-05T11:50:18 ['/Users/dnanto/GitHub/ffbio/ffbio/ffdb.py', 'data/315393']
2020-08-05T11:50:18 txid315393[PORGN] AND biomol_genomic[PROP] NOT gbdiv_pat[PROP] NOT gbdiv_syn[PROP] NOT WGS[KYWD] AND 2020/08/05:3000[MDAT]
2020-08-05T11:50:18 count = 0
```

## ffidx.py

The ```ffidx.py``` script creates and/or queries an indexed set of sequence files created via the [index_db](https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html#index_db) method.


```bash
ffidx.py -h
```

```
usage: ffidx.py [-h] [-filenames FILENAMES [FILENAMES ...]] [-dump]
                [-descriptions] [-entry ENTRY [ENTRY ...]]
                [-entry-batch ENTRY_BATCH] [-fi FI] [-fo FO]
                path

retrieve records from an indexed set of sequence files

positional arguments:
  path                  the sequence flat-file or index path

optional arguments:
  -h, --help            show this help message and exit
  -filenames FILENAMES [FILENAMES ...]
                        the list of sequence files to index (default: None)
  -dump                 the flag to dump all of the records (default: False)
  -descriptions, -headers
                        the flag to only output the descriptions (default:
                        False)
  -entry ENTRY [ENTRY ...]
                        the accessions to retrieve (default: None)
  -entry-batch ENTRY_BATCH
                        the file of accessions to retrieve (default: None)
  -fi FI                the sequence file format (input) (default: fasta)
  -fo FO                the sequence file format (output) (default: fasta)
```

### Ex-1

Create an index using multiple multi-GenBank files and query some accessions.


```bash
ls data/oantigen.[1-2].gbk.gz | \
  xargs ffidx.py data/oantigen.db -entry AF390573.1 GU576499.1 -fi gb -filenames | \
  grep -A 2 \>
```

```
>AF390573.1 Vibrio cholerae serogroup O37 O-antigen biosynthesis region, partial sequence
GCCATCCCACTCTGTGGTCGCAGAGCAAGCTCCCTCATGGAAAATAGCGTCAATGGGCCC
GAAATCATCACCGGCCATGATCTGAGCTAGGAAGTCATCTCGATCCATATAGTCGGCGAT
--
>GU576499.1 Vibrio cholerae strain CO845 O-antigen biosynthesis gene locus, partial sequence
AAGGCGTCATGGACCCGAAATCATCACCAGCCATGATCTGAGCTAGGAAGTCATCTCGAT
CCATATAGTCGGCGATCTGTAGGTCAACCAGATTTTTGAACTTACGACCATTTTTCAAAT
```

List the index.


```bash
du -sh data/oantigen.db
```

```
 20K	data/oantigen.db
```

### Ex-2

Same thing, except create the index temporarily in memory.


```bash
ls data/oantigen.[1-2].gbk.gz | \
  xargs ffidx.py ":memory:" -entry AF390573.1 GU576499.1 -fi gb -filenames | \
  grep -A 2 \>
```

```
>AF390573.1 Vibrio cholerae serogroup O37 O-antigen biosynthesis region, partial sequence
GCCATCCCACTCTGTGGTCGCAGAGCAAGCTCCCTCATGGAAAATAGCGTCAATGGGCCC
GAAATCATCACCGGCCATGATCTGAGCTAGGAAGTCATCTCGATCCATATAGTCGGCGAT
--
>GU576499.1 Vibrio cholerae strain CO845 O-antigen biosynthesis gene locus, partial sequence
AAGGCGTCATGGACCCGAAATCATCACCAGCCATGATCTGAGCTAGGAAGTCATCTCGAT
CCATATAGTCGGCGATCTGTAGGTCAACCAGATTTTTGAACTTACGACCATTTTTCAAAT
```

### Ex-3

Dump all sequences as GenBank records.


```bash
ffidx.py data/oantigen.db -dump -fo gb | grep ^DEFINITION
```

```
DEFINITION  Vibrio cholerae genes for O-antigen synthesis, strain MO45, complete
DEFINITION  Vibrio cholerae genes for o-antigen synthesis, strain O22, complete
DEFINITION  Vibrio cholerae serogroup O37 O-antigen biosynthesis region, partial
DEFINITION  Vibrio cholerae strain CO603B O-antigen biosynthesis gene locus,
DEFINITION  Vibrio cholerae strain CO545 O-antigen biosynthesis gene locus,
DEFINITION  Vibrio cholerae strain CO845 O-antigen biosynthesis gene locus,
```

## ffcds.py

The ```ffcds.py``` script extract all CDS records from a GenBank file.


```bash
ffcds.py -h
```

```
usage: ffcds.py [-h] file

extract CDS sequence records from a GenBank file

positional arguments:
  file        the sequence file

optional arguments:
  -h, --help  show this help message and exit
```

### Ex-1

Extract all CDS records.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | ffcds.py - | grep -A 2 \> | head
```

```
>BAA33585.1 AB012956.1|[7865:8810](-) n/a
ATGATCATCGTCACTGGCGGCGCTGGCATGATTGGCAGCAATATTATCAAAGCGCTTAAT
GAGCGCGGTATCACAGACATTTTGGTCGTTGATCATTTGAAAAATGGTCGTAAGTTCAAA
--
>BAA33586.1 AB012956.1|[8923:10441](-) n/a
ATGCATAAACCAACCATTTCTAGTGTAATCGCACTTACCCTGTTAGGCTGCGGCGGTGGA
GAAAGCGGCAATTCGGGCAACACAACACCACCGGTTAAGTACTTTAATGTAAGCTTTTTG
--
>BAA33587.1 AB012956.1|[10521:12714](-) n/a
ATGAAAACTGGCCACCGCCCTCTTTTGAATACTTCGCTGTCTTTACTCGGCGTACTGATA
```

## ffgbktax.py

The ```ffgbktax.py``` script generates a taxonomy mapping from a GenBank file suitable for use with the ```makeblastdb``` command.


```bash
ffgbktax.py -h
```

```
usage: ffgbktax.py [-h] file

create a taxid map suitable for makeblastdb from a GenBank file

positional arguments:
  file        the sequence file

optional arguments:
  -h, --help  show this help message and exit
```

### Ex-1

Generate a taxonomy mapping from a GenBank file.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | ffgbktax.py - 
```

```
AB012956.1 666
AB012957.1 666
AF390573.1 185332
GU576497.1 666
GU576498.1 666
GU576499.1 666
```

## ffuseq.py

The ```ffuseq.py``` script computes the unique set of sequences and an optional tab-separated mapping file.


```bash
ffuseq.py -h
```

```
usage: ffuseq.py [-h] [-map MAP] [-fmt FMT] file

compute the unique set of sequences

positional arguments:
  file        the sequence file

optional arguments:
  -h, --help  show this help message and exit
  -map MAP    the path prefix for the output files (default: None)
  -fmt FMT    the sequence file format (default: fasta)
```

### Ex-1

This example creates a file stream with duplicate records. Output the unique set and a mapping file.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | ffuseq.py - -fmt gb -map data/urec.tsv > data/urec.gbk
```

```grep``` Unique record LOCUS tags.


```bash
grep ^LOCUS data/urec.gbk
```

```
LOCUS       AB012956               46721 bp    DNA     linear   BCT 16-OCT-1999
LOCUS       GU576497               30443 bp    DNA     linear   BCT 25-JUL-2016
LOCUS       AB012957               45993 bp    DNA     linear   BCT 16-OCT-1999
LOCUS       GU576498               19487 bp    DNA     linear   BCT 25-JUL-2016
LOCUS       GU576499               26128 bp    DNA     linear   BCT 01-APR-2011
LOCUS       AF390573               27552 bp    DNA     linear   BCT 08-MAY-2002
```

```cat``` the unique record mapping file.


```bash
cat data/urec.tsv
```

```
idx	key	id	description	length
1	Wp7X+fuoYSmh9S+8KAOlVOjmMJM	AB012956.1	Vibrio cholerae genes for O-antigen synthesis, strain MO45, complete cds	46721
2	aXZq6FaHBW/pvgXXUJkH1mXCBMs	GU576497.1	Vibrio cholerae strain CO603B O-antigen biosynthesis gene locus, partial sequence	30443
3	gDDKQoXOexqsX3QTZMKiTk3PRI0	AB012957.1	Vibrio cholerae genes for o-antigen synthesis, strain O22, complete cds	45993
4	m3SE/AUt/Wg2D2VfZ92MG7otRAs	GU576498.1	Vibrio cholerae strain CO545 O-antigen biosynthesis gene locus, partial sequence	19487
5	oAPc1vtYNp/vJduAtaFEk7xP69E	GU576499.1	Vibrio cholerae strain CO845 O-antigen biosynthesis gene locus, partial sequence	26128
6	sviBFiv+KU0rRAAD1MBHU32gkac	AF390573.1	Vibrio cholerae serogroup O37 O-antigen biosynthesis region, partial sequence	27552
```

