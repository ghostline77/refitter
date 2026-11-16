//! Inverse-skin corrections for constraints measured in representative poses.
use crate::surface::Surface;

/// Converts contact vectors to a smooth outward shell adjustment. This
/// preserves the garment's local shape instead of moving neighbouring vertices
/// in unrelated directions around folds and openings.
#[must_use]
pub fn smooth_normal_offsets(
    normals: &[[f32; 3]],
    contacts: &[Vec3],
    indices: &[u16],
    movable: &[bool],
) -> Vec<f32> {
    let mut offsets: Vec<_> = normals
        .iter()
        .zip(contacts)
        .enumerate()
        .map(|(i, (normal, contact))| {
            if !movable.get(i).copied().unwrap_or(false) {
                0.
            } else {
                contact
                    .dot(Vec3::from_array(*normal).normalize_or_zero())
                    .clamp(0., 0.01)
            }
        })
        .collect();
    for _ in 0..4 {
        let mut sums = vec![0.; offsets.len()];
        let mut counts = vec![0u32; offsets.len()];
        for triangle in indices.chunks_exact(3) {
            for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                let a = usize::from(triangle[a]);
                let b = usize::from(triangle[b]);
                if a >= offsets.len() || b >= offsets.len() {
                    continue;
                }
                sums[a] += offsets[b];
                counts[a] += 1;
                sums[b] += offsets[a];
                counts[b] += 1;
            }
        }
        for i in 0..offsets.len() {
            if !movable.get(i).copied().unwrap_or(false) {
                offsets[i] = 0.;
            } else if counts[i] > 0 {
                offsets[i] = 0.5 * offsets[i] + 0.5 * sums[i] / counts[i] as f32;
            }
        }
    }
    offsets
}
/// Keeps source-coincident vertices together across UV/material splits. A fixed
/// boundary wins over movable copies, preserving native cap continuity.
pub fn cohere_points(
    points: &mut [Vec<[f32; 3]>],
    seams: &[Vec<(usize, usize)>],
    movable: &[Vec<bool>],
) {
    for seam in seams {
        let valid: Vec<_> = seam
            .iter()
            .copied()
            .filter(|&(c, v)| {
                points.get(c).is_some_and(|p| v < p.len())
                    && movable.get(c).is_some_and(|m| v < m.len())
            })
            .collect();
        if valid.len() < 2 {
            continue;
        }
        let fixed: Vec<_> = valid
            .iter()
            .copied()
            .filter(|&(c, v)| !movable[c][v])
            .collect();
        let anchors = if fixed.is_empty() { &valid } else { &fixed };
        let center = anchors
            .iter()
            .map(|&(c, v)| Vec3::from_array(points[c][v]))
            .sum::<Vec3>()
            / anchors.len() as f32;
        for &(c, v) in &valid {
            points[c][v] = center.to_array();
        }
    }
}
use glam::{Mat4, Vec3};

/// Projects a skinned point outside a local receiving surface, returning a bind
/// space correction. Distant/open boundaries are outside the contact envelope.
#[must_use]
pub fn outside_correction(
    point: Vec3,
    bind_to_world: Mat4,
    surface: &Surface,
    clearance: f32,
    reach: f32,
) -> Option<(Vec3, f32)> {
    if bind_to_world.determinant().abs() < 1e-8 {
        return None;
    }
    let (world, depth) = outside_displacement(point, surface, clearance, reach, None)?;
    let correction = bind_to_world.inverse().transform_vector3(world);
    correction.is_finite().then_some((correction, depth))
}

/// World-space contact displacement. A supplied outward normal enables local
/// sheet contacts, rejecting opposing faces and projections beyond open edges.
/// Without that orientation, penetration requires closed-volume evidence.
#[must_use]
pub fn outside_displacement(
    point: Vec3,
    surface: &Surface,
    clearance: f32,
    reach: f32,
    face_normal: Option<Vec3>,
) -> Option<(Vec3, f32)> {
    let (nearest, normal) = surface.closest(point)?;
    let signed = (point - nearest).dot(normal);
    // Preserve native spacing when surfaces are already separated. The envelope
    // repairs intersections, rather than inflating every near-contact region.
    if signed >= -0.001 {
        return None;
    }
    let depth = clearance - signed;
    if depth <= 0.0005
        || point.distance(nearest) > reach
        || !depth.is_finite()
        || face_normal.is_some_and(|face| face.dot(normal) < 0.35)
    {
        return None;
    }
    // Oriented sheets are valid local barriers without enclosing a volume, but
    // points beyond their open edges must not be pulled across the opening.
    if face_normal.is_some() && ((point - nearest) - normal * signed).length() > 0.003 {
        return None;
    }
    // Unoriented points need volume evidence; oriented sheets were checked
    // against matching normals and their open boundaries above.
    if face_normal.is_none() && depth > clearance + 0.0005 && !surface.contains(point) {
        return None;
    }
    Some((normal * depth.min(0.015), depth))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Triangle;
    #[test]
    fn shared_surface_points_keep_the_fixed_cap_closed() {
        let mut points = vec![vec![[0.01, 0., 0.]], vec![[0.02, 0., 0.]]];
        cohere_points(
            &mut points,
            &[vec![(0, 0), (1, 0)]],
            &[vec![false], vec![true]],
        );
        assert_eq!(points[0][0], points[1][0]);
        assert_eq!(points[0][0], [0.01, 0., 0.]);
    }
    #[test]
    fn contact_field_is_outward_smooth_and_keeps_fixed_boundaries() {
        let field = smooth_normal_offsets(
            &[[1., 0., 0.]; 3],
            &[-Vec3::X, Vec3::X * 0.005, Vec3::ZERO],
            &[0, 1, 2],
            &[true, true, false],
        );
        assert!(field.iter().all(|v| *v >= 0. && *v <= 0.005));
        assert!(field[0] > 0. && field[1] < 0.005);
        assert_eq!(field[2], 0.);
    }

    #[test]
    fn a_point_beyond_an_open_edge_is_not_a_volume_penetration() {
        let surface = Surface::new(&[Triangle([
            Vec3::new(0., -1., -1.),
            Vec3::new(0., 1., -1.),
            Vec3::new(0., 0., 1.),
        ])]);
        assert!(
            outside_displacement(Vec3::new(-0.005, 1.01, -1.), &surface, 0.003, 0.04, None)
                .is_none()
        );
    }

    #[test]
    fn oriented_open_sheet_contacts_preserve_edges_and_opposite_faces() {
        let surface = Surface::new(&[Triangle([
            Vec3::new(0., -1., -1.),
            Vec3::new(0., 1., -1.),
            Vec3::new(0., 0., 1.),
        ])]);
        let point = Vec3::new(-0.005, 0., 0.);
        let (delta, _) = outside_displacement(point, &surface, 0.003, 0.04, Some(Vec3::X)).unwrap();
        assert!((point + delta).x >= 0.0029);
        assert!(outside_displacement(point, &surface, 0.003, 0.04, Some(-Vec3::X)).is_none());
        assert!(
            outside_displacement(
                Vec3::new(-0.005, 1.01, -1.),
                &surface,
                0.003,
                0.04,
                Some(Vec3::X)
            )
            .is_none()
        );
    }
    #[test]
    fn inverse_skin_correction_clears_the_posed_surface() {
        let vertices = [
            Vec3::new(-0.2, -1., -1.),
            Vec3::new(0., -1., -1.),
            Vec3::new(0., 1., -1.),
            Vec3::new(-0.2, 1., -1.),
            Vec3::new(-0.2, -1., 1.),
            Vec3::new(0., -1., 1.),
            Vec3::new(0., 1., 1.),
            Vec3::new(-0.2, 1., 1.),
        ];
        let faces = [
            [0, 3, 2],
            [0, 2, 1],
            [4, 5, 6],
            [4, 6, 7],
            [0, 4, 7],
            [0, 7, 3],
            [1, 2, 6],
            [1, 6, 5],
            [0, 1, 5],
            [0, 5, 4],
            [3, 7, 6],
            [3, 6, 2],
        ];
        let surface = Surface::new(&faces.map(|face| Triangle(face.map(|i| vertices[i]))));
        assert!(surface.contains(Vec3::new(-0.005, 0., 0.)));
        assert!(!surface.contains(Vec3::new(0.005, 0., 0.)));
        let transform = Mat4::from_scale_rotation_translation(
            Vec3::new(2., 1., 0.5),
            glam::Quat::from_rotation_z(0.7),
            Vec3::new(-0.005, 0., 0.),
        );
        let original = Vec3::ZERO;
        let point = transform.transform_point3(original);
        let (correction, _) = outside_correction(point, transform, &surface, 0.003, 0.04).unwrap();
        let corrected = transform.transform_point3(original + correction);
        assert!((corrected.x - 0.003).abs() < 1e-6);
        assert!(outside_displacement(point, &surface, 0.003, 0.04, Some(Vec3::Y)).is_none());
    }
}
