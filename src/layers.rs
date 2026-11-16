//! Surface constraints for layered skinned garments in a shared bind space.

use glam::Vec3;

use crate::{FitReport, Triangle};

/// One inward path from an authored support morph or a skeleton guide.
#[derive(Clone, Copy)]
pub struct SupportPath {
    pub direction: Vec3,
    pub maximum: f32,
}

/// Fits vertices under higher-scored surfaces along their individual support
/// paths. Back-facing exits are ignored, so vertices already inside a layer
/// cannot be pulled through the far side. No axis or garment category is assumed.
pub fn fit_under(
    positions: &mut [[f32; 3]],
    paths: &[SupportPath],
    surfaces: &[Triangle],
    clearance: f32,
) -> FitReport {
    let mut report = FitReport::default();
    if positions.len() != paths.len() || !clearance.is_finite() || clearance <= 0. {
        return report;
    }
    for (position, path) in positions.iter_mut().zip(paths) {
        let p = Vec3::from_array(*position);
        let direction = path.direction.normalize_or_zero();
        if !path.maximum.is_finite() || path.maximum <= 0. || direction == Vec3::ZERO {
            continue;
        }
        let mut nearest = path.maximum;
        let mut hit = false;
        for &Triangle([a, b, c]) in surfaces {
            if let Some(t) = entering_hit(p, direction, a, b, c)
                && t <= nearest
            {
                nearest = t;
                hit = true;
            }
        }
        if hit {
            let amount = (nearest + clearance).min(path.maximum);
            *position = (p + direction * amount).to_array();
            report.changed_vertices += 1;
            report.maximum_displacement = report.maximum_displacement.max(amount);
        }
    }
    report
}

/// Fits between the fixed character body and higher-scored clothing. The body
/// is a hard lower bound: support compression cannot push cloth through skin.
/// Layer depth reserves room for intervening garments without changing loose
/// regions that already have sufficient clearance.
pub fn fit_between(
    positions: &mut [[f32; 3]],
    paths: &[SupportPath],
    body: &crate::surface::Surface,
    outer: &crate::surface::Surface,
    body_clearance: f32,
    layer_clearance: f32,
) -> FitReport {
    let mut report = FitReport::default();
    if positions.len() != paths.len() {
        return report;
    }
    for (position, path) in positions.iter_mut().zip(paths) {
        let original = Vec3::from_array(*position);
        let inward = path.direction.normalize_or_zero();
        if inward == Vec3::ZERO || path.maximum <= 0. {
            continue;
        }
        let mut point = original;
        let mut expanded = 0.;
        if let Some((surface, normal)) = body.closest(point) {
            let signed = (point - surface).dot(normal);
            if signed < body_clearance && point.distance(surface) < path.maximum {
                expanded = body_clearance - signed;
                point += normal * expanded;
            }
        }
        let body_limit = body
            .distance(point, inward, false)
            .map_or(path.maximum + expanded, |distance| {
                (distance - body_clearance.min(0.001)).max(0.)
            });
        let limit = (path.maximum + expanded).min(body_limit);
        if let Some(distance) = outer.distance(point, inward, false).filter(|d| *d <= limit) {
            point += inward * (distance + layer_clearance).min(limit);
        }
        let displacement = original.distance(point);
        if displacement > 1e-6 {
            *position = point.to_array();
            report.changed_vertices += 1;
            report.maximum_displacement = report.maximum_displacement.max(displacement);
        }
    }
    report
}

/// Fits complete triangles as well as their vertices. Edge midpoints and face
/// centres catch coarse triangles that cross a layer despite valid endpoints.
pub fn fit_mesh_between(
    positions: &mut [[f32; 3]],
    indices: &[u16],
    paths: &[SupportPath],
    body: &crate::surface::Surface,
    outer: &crate::surface::Surface,
    clearance: [f32; 2],
) -> FitReport {
    let original = positions.to_vec();
    if positions.len() != paths.len() {
        return FitReport::default();
    }
    for _ in 0..4 {
        fit_between(positions, paths, body, outer, clearance[0], clearance[1]);
        let mut corrections = vec![Vec3::ZERO; positions.len()];
        for triangle in indices.chunks_exact(3) {
            let ids: [usize; 3] = std::array::from_fn(|i| usize::from(triangle[i]));
            if ids.iter().any(|&i| i >= positions.len()) {
                continue;
            }
            for weights in [[1. / 3.; 3], [0.5, 0.5, 0.], [0., 0.5, 0.5], [0.5, 0., 0.5]] {
                let mut point = Vec3::ZERO;
                let mut direction = Vec3::ZERO;
                let mut maximum = 0.;
                for i in 0..3 {
                    point += Vec3::from_array(positions[ids[i]]) * weights[i];
                    direction += paths[ids[i]].direction.normalize_or_zero() * weights[i];
                    maximum += paths[ids[i]].maximum * weights[i];
                }
                // Native closing caps can legitimately span the body interior.
                // Interior samples constrain clothing overlap only; skin
                // clearance is enforced on the actual garment vertices.
                let direction = direction.normalize_or_zero();
                let Some(distance) = outer
                    .distance(point, direction, false)
                    .filter(|d| *d <= maximum)
                else {
                    continue;
                };
                let limit = body
                    .distance(point, direction, false)
                    .map_or(maximum, |d| (d - 0.001).max(0.));
                let correction = direction * (distance + clearance[1]).min(maximum).min(limit);
                let divisor = weights.iter().map(|w| w * w).sum::<f32>();
                for i in 0..3 {
                    let delta = correction * (weights[i] / divisor);
                    if delta.length_squared() > corrections[ids[i]].length_squared() {
                        corrections[ids[i]] = delta;
                    }
                }
            }
        }
        let mut moved = false;
        for (point, delta) in positions.iter_mut().zip(corrections) {
            if delta.length_squared() > 1e-10 {
                *point = (Vec3::from_array(*point) + delta).to_array();
                moved = true;
            }
        }
        for ((point, original), path) in positions.iter_mut().zip(&original).zip(paths) {
            let origin = Vec3::from_array(*original);
            let delta = Vec3::from_array(*point) - origin;
            *point = (origin + delta.clamp_length_max(path.maximum)).to_array();
        }
        if !moved {
            break;
        }
    }
    fit_between(positions, paths, body, outer, clearance[0], clearance[1]);
    let mut report = FitReport::default();
    for (old, point) in original.iter().zip(positions) {
        let distance = Vec3::from_array(*old).distance(Vec3::from_array(*point));
        if distance > 1e-6 {
            report.changed_vertices += 1;
            report.maximum_displacement = report.maximum_displacement.max(distance);
        }
    }
    report
}

/// Distance to the first entering surface, for native coverage decisions.
pub fn entering_distance(origin: Vec3, direction: Vec3, surfaces: &[Triangle]) -> Option<f32> {
    surfaces
        .iter()
        .filter_map(|&Triangle([a, b, c])| entering_hit(origin, direction, a, b, c))
        .min_by(f32::total_cmp)
}

/// Checks a closed receiving volume along three independent directions.
/// Coincident triangle-edge hits count once. Signed crossings allow overlapping
/// closed layers (for example, a shoe upper and sole) without cancelling their
/// coverage as an even/odd parity test would.
pub fn encloses(point: Vec3, surfaces: &[Triangle]) -> bool {
    [Vec3::X, Vec3::Y, Vec3::Z].into_iter().all(|direction| {
        let mut hits: Vec<_> = surfaces
            .iter()
            .filter_map(|&Triangle([a, b, c])| {
                entering_hit(point, direction, a, b, c)
                    .map(|distance| (distance, -1))
                    .or_else(|| {
                        entering_hit(point, direction, c, b, a).map(|distance| (distance, 1))
                    })
            })
            .filter(|(distance, _)| *distance > 1e-5)
            .collect();
        hits.sort_by(|a, b| a.0.total_cmp(&b.0));
        hits.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-5 && a.1 == b.1);
        hits.iter().map(|(_, crossing)| crossing).sum::<i32>() > 0
    })
}

pub(crate) fn entering_hit(p: Vec3, d: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Option<f32> {
    let e1 = b - a;
    let e2 = c - a;
    let h = d.cross(e2);
    let determinant = e1.dot(h);
    // Positive determinant is a ray entering the outward-facing surface.
    if determinant <= 1e-9 {
        return None;
    }
    let inverse = determinant.recip();
    let s = p - a;
    let u = s.dot(h) * inverse;
    if !(-1e-5..=1.00001).contains(&u) {
        return None;
    }
    let q = s.cross(e1);
    let v = d.dot(q) * inverse;
    if v < -1e-5 || u + v > 1.00001 {
        return None;
    }
    let t = e2.dot(q) * inverse;
    (t >= 0. && t.is_finite()).then_some(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn triangle_interior_contact_is_resolved_when_its_vertices_miss_the_surface() {
        let plane = |size: f32, x| {
            Triangle([
                Vec3::new(x, -size, -size),
                Vec3::new(x, size, -size),
                Vec3::new(x, 0., size),
            ])
        };
        let body = crate::surface::Surface::new(&[plane(10., 0.)]);
        let outer = crate::surface::Surface::new(&[plane(0.2, 0.02)]);
        let mut points = [[0.03, -2., -2.], [0.03, 2., -2.], [0.03, 0., 4.]];
        let paths = [SupportPath {
            direction: -Vec3::X,
            maximum: 0.1,
        }; 3];
        let report = fit_mesh_between(
            &mut points,
            &[0, 1, 2],
            &paths,
            &body,
            &outer,
            [0.002, 0.006],
        );
        assert_eq!(report.changed_vertices, 3);
        assert!(points.iter().map(|p| p[0]).sum::<f32>() / 3. < 0.02);
        assert!(points.iter().all(|p| p[0] > 0.));
    }
    #[test]
    fn layer_compression_never_crosses_the_fixed_body_surface() {
        let plane = |x| {
            Triangle([
                Vec3::new(x, -1., -1.),
                Vec3::new(x, 1., -1.),
                Vec3::new(x, 0., 1.),
            ])
        };
        let body = crate::surface::Surface::new(&[plane(0.)]);
        let outer = crate::surface::Surface::new(&[plane(0.02)]);
        let mut points = [[0.03, 0., 0.], [-0.01, 0., 0.]];
        let paths = [SupportPath {
            direction: -Vec3::X,
            maximum: 0.1,
        }; 2];
        fit_between(&mut points, &paths, &body, &outer, 0.004, 0.002);
        assert!((points[0][0] - 0.018).abs() < 1e-5);
        assert!((points[1][0] - 0.004).abs() < 1e-5);
    }
    #[test]
    fn support_paths_work_in_arbitrary_axes_without_pulling_inside_vertices_out() {
        let surfaces = [Triangle([
            Vec3::new(0., -1., -1.),
            Vec3::new(0., 1., -1.),
            Vec3::new(0., 0., 1.),
        ])];
        let mut points = [[0.02, 0., 0.], [-0.02, 0., 0.], [0.02, 3., 0.]];
        let paths = [SupportPath {
            direction: -Vec3::X,
            maximum: 0.1,
        }; 3];
        let report = fit_under(&mut points, &paths, &surfaces, 0.003);
        assert_eq!(report.changed_vertices, 1);
        assert!((points[0][0] + 0.003).abs() < 1e-6);
        assert_eq!(points[1], [-0.02, 0., 0.]);
        assert_eq!(points[2], [0.02, 3., 0.]);
    }
}
