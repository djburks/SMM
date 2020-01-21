# SMM
Single-Order Markov Model (SMM) Builder

A small application for generating single-order markov models for metagenomic classification.  The base program accepts metagenomic read files (in FASTA format), a text directory of genomic FASTA files to model, and user-defined order between 2 and 12.  The log probability of each read originating from each supplied genome is printed to stout (or redirected with the supplied Python script).  

Metagenomic reads are fully loaded into memory, and converted to markov-chain indices for rapid calculation of probabilities. Markov chains are built directly from genomic FASTA files, and kept in memory until all reads have been classified.  There is no database construction step.

Models are pseudo-count supplemented.  An initial count of 4 is assumed for all kmers (user-defined order k), representing each of the four possible combinations of a given kmer and nucleotide {A,T,G,C} at position k+1.

## Quick Start

The program can quickly be compiled with GCC.  For maximum speed, the -Ofast parameter can be specified during compilation without any effect on calculated probabilities.
```
g++ -Ofast smm.cpp -o SMM
```
The supplied Python script is an optional intermediary that generates a user-defined output file with header.  This script can also build the required genomic FASTA directory file by setting the g_dir parameter to the folder containing the genomic FASTA files you wish to model.  If the smm program is in a different directory than the script, please update the smm parameter of smm.py to this location.

Alternatively, the program can be invoked directly.  A directory listing of genomes to be modeled must be supplied.  The full path of each genome should be provided on each line of the text file without spaces.  An example Genomes.txt directory has been provided.

### Running with SMM.py
```
python3 SMM.py Metagenomic_Reads.fasta Metagenomic_Read_Scores.tsv 12
```
### Running Directly
```
./smm Genomes.txt Metagenomic_Reads.fasta 12 >> Metagenomic_Read_Scores.tsv
```

