# NEROH

**Ne Estimation by Runs Of Homozygosity**

NEROH estimates effective population size (Ne) over time from genomic data
using the distribution of runs of homozygosity (ROH). It reads PLINK-format
`.ped`/`.map` files, computes ROH statistics, and fits Ne trajectories using
an evolutionary optimization algorithm.

Authors: Enrique Santiago, Armando Caballero and Carlos Köpke

## System Requirements

- A 64-bit Linux, macOS, or Windows (with WSL/MinGW) system.
- Multi-core CPU recommended. NEROH uses OpenMP to parallelize computation
  across chromosomes and individuals, and benefits significantly from
  additional cores.
- Memory usage depends on dataset size. The program uses compile-time array
  limits of up to 1,000,000 loci and 1,000 individuals.

## Software Requirements

- **C++ compiler** with C++11 support and OpenMP (e.g., GCC 4.9+, Clang with
  libomp).
- **GNU Make** (optional, for using the provided Makefile).

On Debian/Ubuntu:

```bash
sudo apt install g++ make
```

On macOS (with Homebrew):

```bash
brew install gcc
```

Note: Apple Clang does not ship with OpenMP support. Use Homebrew GCC
(`g++-13` or later) and update the `CXX` variable in the Makefile accordingly.

## Installation

Clone or download this repository, then build:

```bash
make
```

This produces the `neroh` executable in the project directory.

To compile manually without Make:

```bash
g++ -O2 -std=c++11 -fopenmp neroh_v2.cpp lib/progress.cpp -o neroh
```

To remove build artifacts:

```bash
make clean
```

## Usage
```
neroh [OPTIONS] <file_name>
```

`file_name` is specified **without** the `.ped` or `.map` extension. The
program expects both `<file_name>.ped` and `<file_name>.map` to exist.

### Examples

```bash
# Phased diploid data (default) with MAF threshold 0.05
./neroh -o 2 -M 0.05 mydata

# Custom distance bounds and marker density cutoff
./neroh -l 0.5 -u 20 -c 50 mydata

# Use 8 threads
./neroh -t 8 mydata
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-h` | Print help message | |
| `-t N` | Number of threads | All available |
| `-o N` | Ploidy mode: 0 = unphased diploid, 1 = haploid, 2 = phased diploid | 2 |
| `-X` | X chromosome mode (compresses distances by 2/3) | Off |
| `-m F` | Mutation rate per Morgan | 0 |
| `-e F` | Genotyping error rate per site | 0 |
| `-l F` | Lower bound distance in cM | 1.0 |
| `-u F` | Upper bound distance in cM | 30 |
| `-r F` | Constant recombination rate in cM/Mb (overrides map cM values) | Use map |
| `-b F` | Bin size in cM for ROH estimation | 0.125 |
| `-d F` | Minimum distance between adjacent markers in cM | 0 |
| `-c N` | Maximum number of markers per cM | All |
| `-M F` | Minor allele frequency (MAF) threshold (forced to 0 when `-m` is given) | 0.05 |
| `-a` | Perform only heterozygosity and marker density analysis (no Ne estimation) | Off |
| `-p` | Print analysis to stdout instead of files | Off |

### Input Format

NEROH uses standard PLINK `.map` and `.ped` files:

- **`.map` file**: Four columns per line <E2><80><94> chromosome, SNP name, position in
  cM, position in bp.
- **`.ped` file**: Six metadata columns (family ID, individual ID, paternal
  ID, maternal ID, sex, phenotype) followed by two allele columns per marker.

### Output Files

Given an input file named `mydata`, NEROH produces:

| File | Description |
|------|-------------|
| `mydata_NeROH_SITES.txt` | Per-block heterozygosity and marker density statistics |
| `mydata_NeROH_STATS.txt` | Summary statistics of the analysis |
| `mydata_NeROH_INPUT.txt` | Preprocessed ROH input data |
| `mydata_NeROH_Ne.txt` | Estimated effective population size (Ne) per generation |
| `mydata_NeROH_DISTRIB.txt` | Observed and predicted ROH length distributions |
| `mydata_NEROH_progress.tmp` | Progress report (can be monitored with `cat` during execution) |
