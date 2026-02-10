# Ilesta

Ilesta is a de novo genome assembler for long reads using the Overlap–Layout–Consensus (OLC) approach. It processes long-read sequencing data to detect overlaps, construct an overlap graph, and generate assembly unitigs.

## Installation

Clone and build using Rust:

```
git clone https://github.com/yvlaere/Ilesta.git
cd Ilesta
cargo build --release
```

## Quick Start

To run the full assembly workflow end to end:

```
Usage: Ilesta assemble [OPTIONS] --reads-fq <READS_FQ>

Options:
  -f, --input-paf <INPUT_PAF>
          Input PAF file (optional if --overlaps is provided)
  -l, --min-overlap-length <MIN_OVERLAP_LENGTH>
          Minimum overlap length [default: 2000]
  -c, --min-overlap-count <MIN_OVERLAP_COUNT>
          Minimum overlap count [default: 3]
  -i, --min-percent-identity <MIN_PERCENT_IDENTITY>
          Minimum percent identity [default: 5]
      --overhang-ratio <OVERHANG_RATIO>
          Overhang ratio [default: 0.8]
      --overlaps <OVERLAPS>
          Pre-computed overlaps binary file (optional, if provided skips alignment filtering)
  -r, --reads-fq <READS_FQ>
          Input reads in FASTQ format
  -p, --output-prefix <OUTPUT_PREFIX>
          Output prefix [default: unitigs]
  -o, --output-dir <OUTPUT_DIR>
          Output directory [default: .]
  -h, --help
          Print help
```

## Development

Ilesta is under active development.