# Ilesta

Ilesta is a de novo genome assembler for long reads using the Overlap–Layout–Consensus (OLC) approach. It processes long-read sequencing data to detect overlaps, construct an overlap graph, and generate assembly unitigs.

## Installation

Clone and build using Rust:

git clone https://github.com/yvlaere/Ilesta.git
cd Ilesta
cargo build --release

## Quick Start

To run the full assembly workflow end to end:

Ilesta assemble
  --input-paf raw_overlaps.paf 
  --min-overlap-length
  --min-overlap-count
  --min-percent-identity
  --overhang-ratio
  --overlap-binary
  --reads-fq
  --output-prefix
  --output-dir

## Development

Ilesta is under active development.