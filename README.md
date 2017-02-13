# bioinfo-scripts
Collection of bioinformatics utility scripts, mostly written in Bioperl

## Naming convention
* `fa-*` - scripts for FASTA files
* `fq-*` - scripts for FASTQ files (4 line format only)
* `gb-*` or `genbank*` - scripts for GENBANK files (sometimes EMBL via `--format` option)
* `fx-*` - scripts that work on FASTA or FASTQ
* `blast*` - scripts that work on standard BLAST output
* `gff*` - scripts for GFF format (usually version 3.0 but check `--gffver` option)
* `faa-*` - scripts for proteins / amino-acid FASTA files

## Other scripts
* `gene-puller.pl` - extract and align amplicons from sets of contigs

## Advice
Some of these scripts were written as early as 2004. Use them to get out of
a bind, for educational reasons, or as a basis for something more robust.
