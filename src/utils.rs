/// General functions used across the project

/// Get the reverse-complement of a node (flip trailing '+' <-> '-').
pub fn rc_node(id: &str) -> String {
    if let Some(last) = id.chars().last() {
        if last == '+' {
            let base = &id[..id.len()-1];
            return format!("{}-", base);
        } else if last == '-' {
            let base = &id[..id.len()-1];
            return format!("{}+", base);
        }
    }
    id.to_string()
}