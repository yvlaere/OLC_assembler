pub struct AlignmentFilteringConfig {
    pub input_paf: String,
    pub min_overlap_length: u32,
    pub min_overlap_count: u32,
    pub min_percent_identity: f32,
    pub overhang_ratio: f32,
}

pub struct CreateOverlapGraphConfig {
    pub overlap_binary: String,
}

pub struct SimplifyOverlapGraphConfig {
    pub max_bubble_length: u32,
    pub tip_length: u32,
}

pub struct UnitigConfig {
    pub overlap_graph_binary: String,
    pub reads_fq: String,
    pub output_prefix: String,
    pub output_dir: String,
}