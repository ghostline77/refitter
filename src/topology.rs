//! Surface connectivity across native UV and material seams.
use std::collections::BTreeMap;

/// Groups faces connected through identical source positions. Vertex and index
/// streams are never welded or rewritten. Invalid triangle indices return None.
#[must_use]
pub fn connected_faces(
    positions: &[[f32; 3]],
    triangles: &[[usize; 3]],
) -> Option<Vec<Vec<usize>>> {
    if positions.iter().flatten().any(|value| !value.is_finite())
        || triangles
            .iter()
            .flatten()
            .any(|&index| index >= positions.len())
    {
        return None;
    }
    let mut parents: Vec<_> = (0..positions.len()).collect();
    let mut welded = BTreeMap::new();
    for (i, point) in positions.iter().enumerate() {
        let key = point.map(|value| if value == 0. { 0 } else { value.to_bits() });
        if let Some(&other) = welded.get(&key) {
            join(&mut parents, i, other);
        } else {
            welded.insert(key, i);
        }
    }
    let valid = |triangle: &[usize; 3]| {
        let [a, b, c] = triangle.map(|i| glam::Vec3::from_array(positions[i]));
        (b - a).cross(c - a).length_squared() > 1e-14
    };
    for &[a, b, c] in triangles.iter().filter(|triangle| valid(triangle)) {
        join(&mut parents, a, b);
        join(&mut parents, a, c);
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for (i, triangle) in triangles.iter().enumerate() {
        if !valid(triangle) {
            continue;
        }
        groups
            .entry(root(&mut parents, triangle[0]))
            .or_default()
            .push(i);
    }
    Some(groups.into_values().collect())
}
fn root(parents: &mut [usize], mut vertex: usize) -> usize {
    while parents[vertex] != vertex {
        parents[vertex] = parents[parents[vertex]];
        vertex = parents[vertex];
    }
    vertex
}
fn join(parents: &mut [usize], a: usize, b: usize) {
    let a = root(parents, a);
    let b = root(parents, b);
    if a != b {
        parents[b] = a;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn uv_splits_remain_one_part_and_disjoint_surfaces_remain_separate() {
        let positions = [
            [0., 0., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [1., 0., 0.],
            [0., 1., 0.],
            [1., 1., 0.],
            [0., 3., 0.],
            [1., 3., 0.],
            [0., 4., 0.],
        ];
        let groups = connected_faces(&positions, &[[0, 1, 2], [3, 4, 5], [6, 7, 8]]).unwrap();
        assert_eq!(groups, vec![vec![0, 1], vec![2]]);
        // A masked zero-area face must not connect independent garment parts.
        let groups =
            connected_faces(&positions, &[[0, 1, 2], [3, 4, 5], [6, 7, 8], [2, 6, 6]]).unwrap();
        assert_eq!(groups, vec![vec![0, 1], vec![2]]);
    }
}
