# SMM
Single-Order Markov Model (SMM) Builder

A small application for generating single-order markov models for metagenomic classification.  The base program accepts metagenomic read files (in FASTA format), a text directory of genomic FASTA files to model, and user-defined order.  Please note that high orders (13+) use significant memory.  The log probability of each read originating from each supplied genome is printed to stout (or redirected with the supplied Python script).  Normalized probability scores can be exported for 12th order analyses.  Trying to use an order other than 12 with the norm parameter will result in a warning (if using the supplied Python script), and reversion to raw scores.

Metagenomic reads are fully loaded into memory, and converted to markov-chain indices for rapid calculation of probabilities. Markov chains are built directly from genomic FASTA files, and kept in memory until all reads have been classified.  There is no database construction step.

Models are pseudo-count supplemented.  An initial count of 4 is assumed for all kmers (user-defined order k), representing each of the four possible combinations of a given kmer and nucleotide {A,T,G,C} at position k+1.

## Quick Start

The program can quickly be compiled with GCC.  For maximum speed, the -Ofast parameter can be specified during compilation without any effect on calculated probabilities.
```
g++ -Ofast smm.cpp -o smm
```
The supplied Python script is an optional intermediary that generates a user-defined output file with header.  If the smm program is in a different directory than the script, please update the smm parameter of smm.py to this location (line 8).

Alternatively, the program can be invoked directly.  A directory listing of genomes to be modeled must be supplied.  The full path of each genome should be provided on each line of the text file without spaces.  An example Genomes.txt directory has been provided.  Probability / Scores for each read are generated in order of their appearance in the metagenomic FASTA file.

### Running with SMM.py with Normalization (12th Order)
```
python3 SMM.py Metagenomic_Reads.fasta Metagenomic_Read_Scores.tsv 12 norm
```
### Running Directly with Normalization (12th Order)
```
./smm Genomes.txt Metagenomic_Reads.fasta 12 norm >> Metagenomic_Read_Scores.tsv
```

### Running with SMM.py without Normalization (8th Order)
```
python3 SMM.py Metagenomic_Reads.fasta Metagenomic_Read_Scores.tsv 8 raw
```
### Running Directly without Normalization (8th Order)
```
./smm Genomes.txt Metagenomic_Reads.fasta 8 raw >> Metagenomic_Read_Scores.tsv
```
## A Note on Output

SMM reports the log-probability that a read originated from a genome based on the single-order markov model constructed from all of the information contained within that genome.  There is no magic number or cutoff that dictates whether a score represents an actual taxonomic match.  One genomic model producing a higher score than another simply means that the read being scored is *more likely* to have originated from that model than the other.

Used as an accompaniment to other taxonomic classifiers, such as kmer/alignment-based tools, SMM can help you determine the closest relative within your genomic database from which a read might have originated.  Given the nature of SMM's modeling, the likliehood that a read and its highest scoring model share their lineage increases as you traverse upwards through the taxonomic ranks (Order --> Class --> Phylum), a characteristic not shared by alignment.

It would be disingenuous to provide a binary (yes/no) classification for each read, as the relationship between probability and classification is highly variable dependent on the database being used and the metagenome being classified.  Users are encouraged, however, to test their own databases with simulated fragments to establish a baseline for their conclusions.

For 12th order probability scores, an optional normalization parameter can be invoked by the user.  This divides the probability by the average probability for a read of the same length across all fully sequenced genomes in the NCBI Genbank database.  Values closer to zero indicate a strong relationship between the model and read.  Values above 1 indicate a relationship between the model and read that is worse than that of a randomly generated sequence of nucleotides of the same read length.
