use std::collections::{HashMap, HashSet};

/// Find weakly connected components ("clustered reads") of the graph.
pub fn weakly_connected_components(graph: &crate::create_string_graph::AssemblyGraph,) -> Vec<Vec<String>> {
    
    // Build an undirected adjacency list
    // Because it is undericted, we can move through both incoming and outgoing edges, meaning we can reach all nodes in a component.
    let mut adjacency_list: HashMap<String, Vec<String>> = HashMap::new();

    // ensure every node in graph.nodes has an entry
    for node_id in graph.nodes.keys() {
        adjacency_list.entry(node_id.clone()).or_default();
    }

    // Populate adjacency using outgoing edges (and add reverse edges to make undirected)
    for (source_id, node) in &graph.nodes {
        // Each edge is a tuple of (target_id, length)
        for (target_id, _length) in &node.edges {
            adjacency_list.entry(source_id.clone()).or_default().push(target_id.clone());
            adjacency_list.entry(target_id.clone()).or_default().push(source_id.clone());
        }
    }

    // Find components with BFS/stack
    let mut visited: HashSet<String> = HashSet::new();
    let mut components: Vec<Vec<String>> = Vec::new();

    for start in adjacency_list.keys() {

        // check if already visited
        if visited.contains(start) {
            continue;
        }

        // new component
        let mut component: Vec<String> = Vec::new();
        let mut stack: Vec<String> = vec![start.clone()];
        visited.insert(start.clone());

        while let Some(current) = stack.pop() {
            component.push(current.clone());
            if let Some(neighbors) = adjacency_list.get(&current) {
                for neighbor in neighbors {
                    if !visited.contains(neighbor) {
                        visited.insert(neighbor.clone());
                        stack.push(neighbor.clone());
                    }
                }
            }
        }

        components.push(component);
    }

    components
}

/// Convenience: return component sizes sorted descending
pub fn component_sizes_sorted(graph: &crate::create_string_graph::AssemblyGraph) -> Vec<usize> {
    let mut sizes: Vec<usize> = weakly_connected_components(graph)
        .into_iter()
        .map(|c| c.len())
        .collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    sizes
}