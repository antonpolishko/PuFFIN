# Puffin 

> PuFFIN - A Parameter-free Method to Build Genome-wide Nucleosome Maps from Paired-end Sequencing Data

## Description

PuFFIN is a command line tool for accurate placing of the nucleosomes based on the pair-end reads. It was designed to place non-overlapping nucleosomes using extra length information present in pair-end data-sets.

The tool is written in python. There are no special requirements except for python2.7+.

The software is freely available for academic use. The software is still in development and may contain bugs and not 100% bulletproof.

## Installation

Clone this repo
```
git clone https://github.com/antonpolishko/PuFFIN.git
```

## Usage

### Input

The input file for the tool should contain the reads only for considered chromosome (contig) and is in a BED format obtained by simply parsing the BAM/SAM file (see bam2bedpe.sh for example).

Given that input reads are in input.bam, that contain only reads for particular chromosome

```
./bam2bedpe.sh input.bam > input.bed
python Run.py input.bed
```

The output will be printed in input.bed.nucs using next column format:

```
<Position of the nucleosome center> <width of the peak> <confidence score> <"Fuzziness"> <Level of the curve that was used to detect nucleosome >
```

## License

MIT © [Anton Polishko](http://www.cs.ucr.edu/~polishka)
