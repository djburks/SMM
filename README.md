# WGM
Whole Genome Markov (WGM) Model Builder

A small application and interface for using Whole-Genome Markov (WGM) models for determining the closest relatives of metagenomic reads.  Genomes are converted to Markov models of user-defined order, and used to determine the probability of shared ancestry with the origins of metagenomic reads.

## Quick Start

Copy the wgm binary to your /usr/local/bin or add it to PATH.
Order is an integer ranging from 1 to 12.  Higher orders tend to produce higher accuracy, but are more influenced by pseudocount additions.
```
wgm <genome file> <order> <metagenomic fasta>
```
### File Restrictions and Formatting

By default, wgm outputs scores to stdout in tab-delimited format.
Genome and Metagenomic Fasta files used with wgm must only contain standard nucleotide characters in upper-case {ATGC} in FASTA format.  Read identifiers should not contain tabs within the name.

The included fastaProc.py script can be used to filter out all non-standard characters.

### Using wgm.py

A python script is included that uses wgm to perform bulk analysis of metagenomic reads.  Multiple (filtered) genomes can be included in the Genomes directory.  The script evokes wgm and passes the desired order and output file destination.

```
python3 wgm.py <metagenomic fasta> <output file> <order>
```
Output is a tab-delimited spreadsheet with log10 probability scores for all read (columns) and genome (rows) pairs.
Genomes producing the highest probability for a given read are most likely to share taxonomic origins.

