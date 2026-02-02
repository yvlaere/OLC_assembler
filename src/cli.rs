use clap::{Parser, Subcommand, Args};

#[derive(Parser)]
#[command(name = "Ilesta", version = "0.0", about = "De novo genome assembly for long reads using the OLC approach")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    
    /// Alignment filtering
    AlignmentFiltering(AlignmentFilteringArgs),

    /// Overlap graph construction
    CreateOverlapGraph(CreateOverlapGraphArgs),

    /// Create unitigs
    CreateUnitigs(CreateUnitigsArgs),

    /// Full genome assembly pipeline
    Assemble(AssembleArgs),
}

#[derive(Args)]
pub struct AlignmentFilteringArgs {

    /// Input PAF file
    #[arg(short, long)]
    pub input_paf: String,

    /// Minimum overlap length
    #[arg(short = 'l', long, default_value_t = 2000)]
    pub min_overlap_length: u32,

    /// Minimum overlap count
    #[arg(short = 'c', long, default_value_t = 3)]
    pub min_overlap_count: u32,

    /// Minimum percent identity
    #[arg(short = 'i', long, default_value_t = 5.0)]
    pub min_percent_identity: f32,

    /// Overhang ratio
    #[arg(short = 'h', long, default_value_t = 0.8)]
    pub overhang_ratio: f32,

}

impl From<&AlignmentFilteringArgs> for crate::configs::AlignmentFilteringConfig {
    fn from(args: &AlignmentFilteringArgs) -> Self {
        Self {
            input_paf: args.input_paf,
            min_overlap_length: args.min_overlap_length,
            min_overlap_count: args.min_overlap_count,
            min_percent_identity: args.min_percent_identity,
            overhang_ratio: args.overhang_ratio,
        }
    }
}

#[derive(Args)]
pub struct CreateOverlapGraphArgs {

    /// Overlap binary file
    #[arg(short = 'b', long, default_value = "overlaps.bin")]
    pub overlap_binary: String,
}

impl From<&CreateOverlapGraphArgs> for crate::configs::CreateOverlapGraphConfig {
    fn from(args: &CreateOverlapGraphArgs) -> Self {
        Self {
            overlap_binary: args.overlap_binary.clone(),
        }
    }
}

#[derive(Args)]
pub struct CreateUnitigsArgs {
    

#[derive(Args)]
pub struct AssembleArgs {

    /// Input PAF file
    #[arg(short = 'a', long = "alignments")]
    pub input_paf: String,

    /// Input reads in FASTQ format
    #[arg(short = 'r', long = "reads")]
    pub reads_fq: String,

    /// Output prefix
    #[arg(short = 'p', long = "prefix", default_value = "assembly")]
    pub output_prefix: String,

    /// Output directory
    #[arg(short = 'o', long = "out_dir", default_value = ".")]
    pub output_dir: String,
}