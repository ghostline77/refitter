//! Proper noncoplanar triangle crossings for fitted garment diagnostics.
use crate::Triangle;
use glam::Vec3;

fn segment_hit(start: Vec3, end: Vec3, Triangle([a, b, c]): Triangle) -> Option<Vec3> {
    let direction = end - start;
    let ab = b - a;
    let ac = c - a;
    let p = direction.cross(ac);
    let determinant = ab.dot(p);
    if determinant.abs() < 1e-10 {
        return None;
    }
    let delta = start - a;
    let u = delta.dot(p) / determinant;
    let q = delta.cross(ab);
    let v = direction.dot(q) / determinant;
    let t = ac.dot(q) / determinant;
    (u >= -1e-5 && v >= -1e-5 && u + v <= 1.00001 && (-1e-5..=1.00001).contains(&t))
        .then_some(start + direction * t)
}
/// Returns a proper intersection segment. Coplanar overlaps and isolated point
/// touches are excluded, as are invalid and degenerate input triangles.
#[must_use]
pub fn segment(a: Triangle, b: Triangle) -> Option<[Vec3; 2]> {
    if a.0.iter().chain(&b.0).any(|point| !point.is_finite()) {
        return None;
    }
    let mut points = Vec::new();
    for (source, target) in [(a, b), (b, a)] {
        for (i, j) in [(0, 1), (1, 2), (2, 0)] {
            if let Some(point) = segment_hit(source.0[i], source.0[j], target)
                && !points
                    .iter()
                    .any(|old: &Vec3| old.distance_squared(point) < 1e-10)
            {
                points.push(point);
            }
        }
    }
    let mut longest = None;
    let mut length = 1e-8;
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            let d = points[i].distance_squared(points[j]);
            if d > length {
                length = d;
                longest = Some([points[i], points[j]]);
            }
        }
    }
    longest
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn proper_crossing_is_distinct_from_disjoint_and_coplanar_surfaces() {
        let a = Triangle([
            Vec3::new(-1., -1., 0.),
            Vec3::new(1., -1., 0.),
            Vec3::new(0., 1., 0.),
        ]);
        let b = Triangle([
            Vec3::new(0., -0.5, -1.),
            Vec3::new(0., -0.5, 1.),
            Vec3::new(0., 0.5, 0.),
        ]);
        let line = segment(a, b).unwrap();
        assert!((line[0].distance(line[1]) - 1.).abs() < 1e-5);
        assert!(segment(a, Triangle(b.0.map(|p| p + Vec3::X * 3.))).is_none());
        assert!(segment(a, a).is_none());
    }
}
