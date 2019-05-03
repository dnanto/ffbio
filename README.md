---
title: "README"
author: "dnanto"
date: "May 03, 2019"
output: 
  html_document: 
    keep_md: yes
    toc: yes
    toc_depth: 2
---



# Abstract

**ffdb** is a collection of scripts to work with flat-file sequence databases and help me complete my thesis.

Note: this repository is experimental...

# Dependencies

* The Pipfile lists the Python requirements.
* [edirect](ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/) should be in PATH.
* For convenience, add the script root directory to PATH.
* Also, it is convenient to know the NCBI [search field descriptions](https://www.ncbi.nlm.nih.gov/books/NBK49540/).

# Scripts

## ffdb.py

The ```ffdb.py``` script searches NCBI and downloads the entire result set. The destination folder stores the [BGZip](http://www.htslib.org/doc/bgzip.html)-compressed flat-files. The [index_db](https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html#index_db) method generates the [SQLite](https://www.sqlite.org/index.html) database in the same directory. The script also modifies the meta_data table to preserve NCBI search parameters.


```bash
ffdb.py -h
```

```
usage: ffdb.py [-h] [-fmt FMT] [-db DB] [-term TERM] [-no-mdat] [-email EMAIL]
               [-post-size POST_SIZE] [-cache CACHE [CACHE ...]]
               [-cache-fmt CACHE_FMT] [-redo]
               repo

update a repository of indexed sequence files

positional arguments:
  repo                  the target file to create/update

optional arguments:
  -h, --help            show this help message and exit
  -fmt FMT, --fmt FMT, -format FMT, --format FMT
                        the sequence file format (default: fasta)
  -db DB, --db DB, -database DB, --database DB
                        the NCBI database (default: nuccore)
  -term TERM, --term TERM
                        the NCBI query term (default: None)
  -no-mdat, --no-mdat   flag to disable adding last modified date to the query
                        (default: False)
  -email EMAIL, --email EMAIL
                        the e-mail to identify yourself to NCBI (for
                        politeness reasons) (default: )
  -post-size POST_SIZE, --post-size POST_SIZE
                        the number of records to post at a time (default:
                        1000)
  -cache CACHE [CACHE ...], --cache CACHE [CACHE ...]
                        the cache of sequence files to search prior to
                        querying NCBI (default: None)
  -cache-fmt CACHE_FMT, --cache-fmt CACHE_FMT
                        the sequence file format of the cache (default: fasta)
  -redo, --redo         the flag to delete everything and redo (default:
                        False)
```

### Ex-1

Initialize a repository of indexed GenBank files. Search NCBI for all *Mimivirus* genomic DNA sequences.


```bash
ffdb.py data/315393 -term 'txid315393[PORGN] AND biomol_genomic[PROP]' -fmt gbk -post 100
```

```
term: txid315393[PORGN] AND biomol_genomic[PROP]
count:  246
new:  246
cached: 0
download: 246
download: KF493705.1 MG779323.1 MG779401.1 JQ063128.1 KT321961.1 ...
download: KT595683.1 MG779376.1 MH046814.1 MG779385.1 JN885991.1 ...
download: KT595677.1 KT599914.1 MG779340.1 MG779300.1 MG779391.1 ...
```

List the results in the directory.


```bash
du -sh data/315393.*
```

```
8.0M	data/315393.0.gbk.bgz
 11M	data/315393.1.gbk.bgz
6.0M	data/315393.2.gbk.bgz
 36K	data/315393.idx
```

### Ex-2

Update the same repository.


```bash
ffdb.py data/315393
```

```
term: txid315393[PORGN] AND biomol_genomic[PROP] AND 2019/05/03:3000[MDAT]
count:  0
new:  0
cached: 0
download: 0
```

### Ex-2

Update the same repository, but disable the modified date limiter.


```bash
ffdb.py data/315393 -no-mdat
```

```
term: txid315393[PORGN] AND biomol_genomic[PROP]
count:  246
new:  0
cached: 0
download: 0
```

## ffidx.py

The ```ffidx.py``` script creates and/or queries an indexed set of sequence files created via the [index_db](https://biopython.org/DIST/docs/api/Bio.SeqIO-module.html#index_db) method.


```bash
ffidx.py -h
```

```
usage: ffidx.py [-h] [-filenames FILENAMES [FILENAMES ...]] [-all]
                [-accessions KEYS [KEYS ...]] [-no-version] [-fmt-idx FMT_IDX]
                [-fmt-out FMT_OUT]
                index

retrieve records from an indexed set of sequence files

positional arguments:
  index                 the SQLite index

optional arguments:
  -h, --help            show this help message and exit
  -filenames FILENAMES [FILENAMES ...], --filenames FILENAMES [FILENAMES ...]
                        the list of sequence files to index (default: None)
  -all, -all, -dump, --dump
                        the flag to dump all of the records (default: False)
  -accessions KEYS [KEYS ...], --accessions KEYS [KEYS ...]
                        the accessions to retrieve (default: None)
  -no-version, --no-version
                        the flag to indicate that the accessions are missing a
                        version (default: False)
  -fmt-idx FMT_IDX, --fmt-idx FMT_IDX
                        the sequence file format of the indexed files,
                        optional if reloading (default: None)
  -fmt-out FMT_OUT, --fmt-out FMT_OUT
                        the sequence file format (output) (default: fasta)
```

### Ex-1

Create an index using multiple multi-GenBank files and query some accessions.


```bash
ls data/oantigen.[1-2].gbk.gz | \
  xargs ffidx.py data/oantigen.idx -acc AF390573.1 GU576499.1 -fmt-idx gb -filenames | \
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
du -sh data/oantigen.idx
```

```
 20K	data/oantigen.idx
```

### Ex-2

Same thing, except create the index temporarily in memory.


```bash
ls data/oantigen.[1-2].gbk.gz | \
  xargs ffidx.py ":memory:" -acc AF390573.1 GU576499.1 -fmt-idx gb -filenames | \
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
ffidx.py data/oantigen.idx -dump -fmt-o gb | grep ^DEFINITION
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
>lcl|BAA33585.1 AB012956.1|n/a|[7865:8810](-)
ATGATCATCGTCACTGGCGGCGCTGGCATGATTGGCAGCAATATTATCAAAGCGCTTAAT
GAGCGCGGTATCACAGACATTTTGGTCGTTGATCATTTGAAAAATGGTCGTAAGTTCAAA
--
>lcl|BAA33586.1 AB012956.1|n/a|[8923:10441](-)
ATGCATAAACCAACCATTTCTAGTGTAATCGCACTTACCCTGTTAGGCTGCGGCGGTGGA
GAAAGCGGCAATTCGGGCAACACAACACCACCGGTTAAGTACTTTAATGTAAGCTTTTTG
--
>lcl|BAA33587.1 AB012956.1|n/a|[10521:12714](-)
ATGAAAACTGGCCACCGCCCTCTTTTGAATACTTCGCTGTCTTTACTCGGCGTACTGATA
```

## fffilter.py

The ```fffilter.py``` script filters for records with matching regex patterns and/or length bounds.


```bash
fffilter.py -h
```

```
usage: fffilter.py [-h] [-fmt FMT] [-fmt-o FMT_O] [-pattern PATTERN]
                   [-length LENGTH] [-percentage]
                   file

filter sequence records by header and/or length

positional arguments:
  file                  the sequence file

optional arguments:
  -h, --help            show this help message and exit
  -fmt FMT, --fmt FMT, -format FMT, --format FMT
                        the sequence file format (input) (default: fasta)
  -fmt-o FMT_O, --fmt-o FMT_O
                        the sequence file format (output) (default: fasta)
  -pattern PATTERN      the regex pattern to search headers (default: None)
  -length LENGTH        the sequence length (default: None)
  -percentage, --percentage
                        the flag to use percentage bounds (default: False)
```

### Length Syntax

The length syntax is length:+:-, where length is the sequence length in bp, + is the upper bound, and - is the lower bound. If the ```-percentage``` flag is set then the bounds are interpreted as percentages in the range [0, 100].

### Ex-1

Filter sequences based on a regex pattern.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | fffilter.py - -fmt gb -pattern " CO[0-9]+ " | grep \>
```

```
>GU576498.1 Vibrio cholerae strain CO545 O-antigen biosynthesis gene locus, partial sequence
>GU576499.1 Vibrio cholerae strain CO845 O-antigen biosynthesis gene locus, partial sequence
```

### Ex-2

Filter sequences based on length bounds: 30000 bp +/- 10%.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | fffilter.py - -fmt gb -len 30000:10:10 -per | grep \>
```

```
>AF390573.1 Vibrio cholerae serogroup O37 O-antigen biosynthesis region, partial sequence
>GU576497.1 Vibrio cholerae strain CO603B O-antigen biosynthesis gene locus, partial sequence
```

### Ex-3

Filter sequences based on length bounds: 45000 bp +/- 5000.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | fffilter.py - -fmt gb -len 45000:5000:5000 | grep \>
```

```
>AB012956.1 Vibrio cholerae genes for O-antigen synthesis, strain MO45, complete cds
>AB012957.1 Vibrio cholerae genes for o-antigen synthesis, strain O22, complete cds
```

## fflen.py

The ```fflen.py``` script computes the length for each sequence in the flat-file.


```bash
fflen.py -h
```

```
usage: fflen.py [-h] [-fmt FMT] [-separator SEP] [-summary] file

compute the length of each record in the file

positional arguments:
  file                  the sequence file

optional arguments:
  -h, --help            show this help message and exit
  -fmt FMT, --fmt FMT, -format FMT, --format FMT
                        the sequence file format (default: fasta)
  -separator SEP, --separator SEP
                        the table delimiter, default is the tab character
                        (default: )
  -summary, --summary, -stats, --stats
                        the flag to compute basic summary statistics (default:
                        False)
```

### Ex-1

List the length of each record in the file.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | fflen.py - -fmt gb
```

```
AB012956.1	Vibrio cholerae genes for O-antigen synthesis, strain MO45, complete cds	46721
AB012957.1	Vibrio cholerae genes for o-antigen synthesis, strain O22, complete cds	45993
AF390573.1	Vibrio cholerae serogroup O37 O-antigen biosynthesis region, partial sequence	27552
GU576497.1	Vibrio cholerae strain CO603B O-antigen biosynthesis gene locus, partial sequence	30443
GU576498.1	Vibrio cholerae strain CO545 O-antigen biosynthesis gene locus, partial sequence	19487
GU576499.1	Vibrio cholerae strain CO845 O-antigen biosynthesis gene locus, partial sequence	26128
```

### Ex-2

Compute the length summary statistics.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | fflen.py - -fmt gb -stat
```

```
mean	32720.666666666668
median	28997.5
std	10187.085233547208
min	19487
max	46721
```

## ffgbk2taxmap.py

The ```ffgbk2taxmap.py``` script generates a taxonomy mapping from a GenBank file suitable for use with the ```makeblastdb``` command.


```bash
ffgbk2taxmap.py -h
```

```
usage: ffgbk2taxmap.py [-h] file

create a taxid map suitable for makeblastdb from a GenBank file

positional arguments:
  file        the sequence file

optional arguments:
  -h, --help  show this help message and exit
```

### Ex-1

Generate a taxonomy mapping from a GenBank file.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | ffgbk2taxmap.py - 
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
usage: ffuseq.py [-h] [-fmt FMT] [-map MAP] file

compute the unique set of sequences

positional arguments:
  file                 the sequence file

optional arguments:
  -h, --help           show this help message and exit
  -fmt FMT, --fmt FMT
  -map MAP, --map MAP  the path prefix for the output files (default: None)
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
LOCUS       AB012957               45993 bp    DNA     linear   BCT 16-OCT-1999
LOCUS       AF390573               27552 bp    DNA     linear   BCT 08-MAY-2002
LOCUS       GU576497               30443 bp    DNA     linear   BCT 25-JUL-2016
LOCUS       GU576498               19487 bp    DNA     linear   BCT 25-JUL-2016
LOCUS       GU576499               26128 bp    DNA     linear   BCT 01-APR-2011
```

```cat``` the unique record mapping file.


```bash
cat data/urec.tsv
```

```
idx	key	id	description
1	Wp7X+fuoYSmh9S+8KAOlVOjmMJM	AB012956.1	Vibrio cholerae genes for O-antigen synthesis, strain MO45, complete cds
2	gDDKQoXOexqsX3QTZMKiTk3PRI0	AB012957.1	Vibrio cholerae genes for o-antigen synthesis, strain O22, complete cds
3	sviBFiv+KU0rRAAD1MBHU32gkac	AF390573.1	Vibrio cholerae serogroup O37 O-antigen biosynthesis region, partial sequence
4	aXZq6FaHBW/pvgXXUJkH1mXCBMs	GU576497.1	Vibrio cholerae strain CO603B O-antigen biosynthesis gene locus, partial sequence
5	m3SE/AUt/Wg2D2VfZ92MG7otRAs	GU576498.1	Vibrio cholerae strain CO545 O-antigen biosynthesis gene locus, partial sequence
6	oAPc1vtYNp/vJduAtaFEk7xP69E	GU576499.1	Vibrio cholerae strain CO845 O-antigen biosynthesis gene locus, partial sequence
```

## ffgbksrc.py

The ```ffgbksrc``` command queries GenBank metadata from the **source** feature.


```bash
ffgbksrc.py -h
```

```
usage: ffgbksrc.py [-h] [-fields] [-separator SEP] file keys [keys ...]

output a table of source feature metadata from a GenBank file

positional arguments:
  file                  the sequence file
  keys                  the qualifier keys

optional arguments:
  -h, --help            show this help message and exit
  -fields, -fields, -header, --header
                        the flag to output a header (default: False)
  -separator SEP, --separator SEP
                        the table delimiter, the default is a tab character
                        (default: )
```

### Key Syntax

The key syntax is key:default:join, where key is the key qualifier, default is the default value if the key is missing, and join is a character to join multiple values under the same key. The default and join specifiers are optional.

### Ex-1


```bash
gunzip -c data/oantigen.[1-2].gbk.gz | ffgbksrc.py - strain db_xref:0:,
```

```
AB012956.1	MO45	taxon:666
AB012957.1	O22	taxon:666
AF390573.1	1322-69	taxon:185332
GU576497.1	CO603B	taxon:666
GU576498.1	CO545	taxon:666
GU576499.1	CO845	taxon:666
```

## fffnagff2gbk.py

The ```fffnagff2gbk.py``` script combines FASTA and GFF3 files to create GenBank files with feature annotations.


```bash
fffnagff2gbk.py -h
```

```
usage: fffnagff2gbk.py [-h] -fna FNA [FNA ...] -gff GFF [GFF ...] [-out OUT]

combine FASTA and GFF3 files to create GenBank files

optional arguments:
  -h, --help           show this help message and exit
  -fna FNA [FNA ...]   the sequence files (default: None)
  -gff GFF [GFF ...]   the gff files (default: None)
  -out OUT, --out OUT
```

### Ex-1

...

## ffkmer.py

The ```ffkmer``` script k-merizes sequences into separate records.


```bash
ffkmer.py -h
```

```
usage: ffkmer.py [-h] [-fmt FMT] [-step STEP] file size

calculate the set of k-mers

positional arguments:
  file                  the sequence file
  size                  the k-mer size

optional arguments:
  -h, --help            show this help message and exit
  -fmt FMT, --fmt FMT, -format FMT, --format FMT
                        the sequence file format (default: fasta)
  -step STEP            the next position to scan in the sequence (default: 1)
```

### Ex-1

Calculate all 100-mers.


```bash
gunzip -c data/oantigen.[1-2].gbk.gz 2> /dev/null | ffkmer.py - 100 -fmt gb | head | grep -A 2 \>
```

```
>AB012956.1-0
GATCAATTTTTTCATTGGGCGAATTCCTTAGTAAAACCTAATTCGGTCAAGGTATTTTCGAGAGCTCCACGATTTTGCATTACAATCGCTAAGGCATTTT
>AB012956.1-1
ATCAATTTTTTCATTGGGCGAATTCCTTAGTAAAACCTAATTCGGTCAAGGTATTTTCGAGAGCTCCACGATTTTGCATTACAATCGCTAAGGCATTTTT
>AB012956.1-2
TCAATTTTTTCATTGGGCGAATTCCTTAGTAAAACCTAATTCGGTCAAGGTATTTTCGAGAGCTCCACGATTTTGCATTACAATCGCTAAGGCATTTTTT
>AB012956.1-3
CAATTTTTTCATTGGGCGAATTCCTTAGTAAAACCTAATTCGGTCAAGGTATTTTCGAGAGCTCCACGATTTTGCATTACAATCGCTAAGGCATTTTTTC
>AB012956.1-4
AATTTTTTCATTGGGCGAATTCCTTAGTAAAACCTAATTCGGTCAAGGTATTTTCGAGAGCTCCACGATTTTGCATTACAATCGCTAAGGCATTTTTTCC
```

