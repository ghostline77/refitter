//! Bounding-volume acceleration for garment/body surface constraints.
use crate::Triangle;
use glam::Vec3;
struct Node {
    min: Vec3,
    max: Vec3,
    start: usize,
    end: usize,
    children: Option<(usize, usize)>,
}
/// Reusable triangle surface; construct once for all vertex constraints.
pub struct Surface {
    triangles: Vec<Triangle>,
    nodes: Vec<Node>,
}

/// A proper crossing with a triangle retained in the accelerated surface.
#[derive(Debug)]
pub struct SurfaceCrossing {
    /// Receiving surface triangle.
    pub triangle: Triangle,
    /// Endpoints of the noncoplanar intersection segment, in metres.
    pub segment: [Vec3; 2],
}
impl Surface {
    /// Finds proper triangle crossings through the surface's bounding hierarchy.
    /// Coplanar overlap, point touches and invalid queries do not produce hits.
    #[must_use]
    pub fn crossings_triangle(&self, query: Triangle) -> Vec<SurfaceCrossing> {
        if self.nodes.is_empty() || query.0.iter().any(|p| !p.is_finite()) {
            return Vec::new();
        }
        let min = query
            .0
            .into_iter()
            .fold(Vec3::splat(f32::INFINITY), Vec3::min)
            - Vec3::splat(1e-5);
        let max = query
            .0
            .into_iter()
            .fold(Vec3::splat(f32::NEG_INFINITY), Vec3::max)
            + Vec3::splat(1e-5);
        let mut pending = vec![0];
        let mut result = Vec::new();
        while let Some(index) = pending.pop() {
            let node = &self.nodes[index];
            if !node.min.cmple(max).all() || !node.max.cmpge(min).all() {
                continue;
            }
            if let Some((left, right)) = node.children {
                pending.extend([left, right]);
            } else {
                for &triangle in &self.triangles[node.start..node.end] {
                    if let Some(segment) = crate::intersection::segment(query, triangle) {
                        result.push(SurfaceCrossing { triangle, segment });
                    }
                }
            }
        }
        result
    }
    /// Tests interior coverage using signed crossings along three oblique rays.
    /// Overlapping shells contribute winding rather than cancelling by parity.
    #[must_use]
    pub fn contains(&self, point: Vec3) -> bool {
        if self.nodes.is_empty() {
            return false;
        }
        [
            Vec3::new(1., 0.317, 0.173),
            Vec3::new(0.193, 1., 0.417),
            Vec3::new(0.371, 0.137, 1.),
        ]
        .into_iter()
        .filter(|direction| {
            let mut hits = Vec::new();
            self.crossings(0, point, *direction, &mut hits);
            hits.sort_by(|a, b| a.0.total_cmp(&b.0));
            hits.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-5 && a.1 == b.1);
            hits.iter().map(|(_, side)| side).sum::<i32>() > 0
        })
        .count()
            >= 2
    }
    fn crossings(&self, index: usize, p: Vec3, d: Vec3, hits: &mut Vec<(f32, i32)>) {
        let node = &self.nodes[index];
        let mut low = 0_f32;
        let mut high = f32::INFINITY;
        for axis in 0..3 {
            let a = (node.min[axis] - p[axis]) / d[axis];
            let b = (node.max[axis] - p[axis]) / d[axis];
            low = low.max(a.min(b));
            high = high.min(a.max(b));
            if low > high {
                return;
            }
        }
        if let Some((left, right)) = node.children {
            self.crossings(left, p, d, hits);
            self.crossings(right, p, d, hits);
        } else {
            for &Triangle([a, b, c]) in &self.triangles[node.start..node.end] {
                if let Some(distance) =
                    crate::layers::entering_hit(p, d, a, b, c).filter(|t| *t > 1e-5)
                {
                    hits.push((distance, -1));
                }
                if let Some(distance) =
                    crate::layers::entering_hit(p, d, c, b, a).filter(|t| *t > 1e-5)
                {
                    hits.push((distance, 1));
                }
            }
        }
    }
    /// Nearest point and outward normal, with distance-based BVH pruning.
    #[must_use]
    pub fn closest(&self, point: Vec3) -> Option<(Vec3, Vec3)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best = f32::INFINITY;
        let mut result = None;
        self.closest_visit(0, point, &mut best, &mut result);
        result.map(|(point, normal, _)| (point, normal))
    }
    /// Nearest source triangle and point, for transferring surface deformation.
    #[must_use]
    pub fn closest_triangle(&self, point: Vec3) -> Option<(Vec3, Triangle)> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut best = f32::INFINITY;
        let mut result = None;
        self.closest_visit(0, point, &mut best, &mut result);
        result.map(|(point, _, triangle)| (point, triangle))
    }
    fn closest_visit(
        &self,
        index: usize,
        point: Vec3,
        best: &mut f32,
        result: &mut Option<(Vec3, Vec3, Triangle)>,
    ) {
        let node = &self.nodes[index];
        if point.distance_squared(point.clamp(node.min, node.max)) > *best {
            return;
        }
        if let Some((left, right)) = node.children {
            let distance = |index: usize| {
                let n = &self.nodes[index];
                point.distance_squared(point.clamp(n.min, n.max))
            };
            let (first, second) = if distance(left) <= distance(right) {
                (left, right)
            } else {
                (right, left)
            };
            self.closest_visit(first, point, best, result);
            self.closest_visit(second, point, best, result);
        } else {
            for &Triangle([a, b, c]) in &self.triangles[node.start..node.end] {
                let ab = b - a;
                let ac = c - a;
                let normal = ab.cross(ac).normalize_or_zero();
                if normal == Vec3::ZERO {
                    continue;
                }
                let projected = point - normal * (point - a).dot(normal);
                let inside = [(a, b), (b, c), (c, a)]
                    .iter()
                    .all(|&(u, v)| (v - u).cross(projected - u).dot(normal) >= -1e-9);
                let closest = if inside {
                    projected
                } else {
                    [(a, b), (b, c), (c, a)]
                        .into_iter()
                        .map(|(u, v)| {
                            let edge = v - u;
                            u + edge
                                * ((point - u).dot(edge) / edge.length_squared().max(1e-20))
                                    .clamp(0., 1.)
                        })
                        .min_by(|u, v| {
                            point
                                .distance_squared(*u)
                                .total_cmp(&point.distance_squared(*v))
                        })
                        .unwrap_or(a)
                };
                let distance = point.distance_squared(closest);
                if distance < *best {
                    *best = distance;
                    *result = Some((closest, normal, Triangle([a, b, c])));
                }
            }
        }
    }
    #[must_use]
    pub fn new(triangles: &[Triangle]) -> Self {
        let mut result = Self {
            triangles: triangles.to_vec(),
            nodes: Vec::new(),
        };
        if !result.triangles.is_empty() {
            result.build(0, result.triangles.len());
        }
        result
    }
    fn build(&mut self, start: usize, end: usize) -> usize {
        let mut min = Vec3::splat(f32::INFINITY);
        let mut max = Vec3::splat(f32::NEG_INFINITY);
        for Triangle(points) in &self.triangles[start..end] {
            for point in points {
                min = min.min(*point);
                max = max.max(*point);
            }
        }
        let index = self.nodes.len();
        self.nodes.push(Node {
            min,
            max,
            start,
            end,
            children: None,
        });
        if end - start > 12 {
            let span = max - min;
            let axis = if span.x >= span.y && span.x >= span.z {
                0
            } else if span.y >= span.z {
                1
            } else {
                2
            };
            self.triangles[start..end].sort_unstable_by(|a, b| {
                let center = |t: &Triangle| t.0.iter().map(|p| p[axis]).sum::<f32>();
                center(a).total_cmp(&center(b))
            });
            let middle = (start + end) / 2;
            let left = self.build(start, middle);
            let right = self.build(middle, end);
            self.nodes[index].children = Some((left, right));
        }
        index
    }
    /// Nearest outward-surface entry, or exit when `exit` is true, along a unit ray.
    #[must_use]
    pub fn distance(&self, origin: Vec3, direction: Vec3, exit: bool) -> Option<f32> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut nearest = f32::INFINITY;
        self.visit(0, origin, direction, exit, &mut nearest);
        nearest.is_finite().then_some(nearest)
    }
    fn visit(&self, index: usize, p: Vec3, d: Vec3, exit: bool, nearest: &mut f32) {
        let node = &self.nodes[index];
        let mut low = 0_f32;
        let mut high = *nearest;
        for axis in 0..3 {
            if d[axis].abs() < 1e-9 {
                if p[axis] < node.min[axis] || p[axis] > node.max[axis] {
                    return;
                }
                continue;
            }
            let a = (node.min[axis] - p[axis]) / d[axis];
            let b = (node.max[axis] - p[axis]) / d[axis];
            low = low.max(a.min(b));
            high = high.min(a.max(b));
            if low > high {
                return;
            }
        }
        if let Some((left, right)) = node.children {
            self.visit(left, p, d, exit, nearest);
            self.visit(right, p, d, exit, nearest);
        } else {
            for &Triangle([a, b, c]) in &self.triangles[node.start..node.end] {
                let hit = if exit {
                    crate::layers::entering_hit(p, d, c, b, a)
                } else {
                    crate::layers::entering_hit(p, d, a, b, c)
                };
                if let Some(distance) = hit {
                    *nearest = nearest.min(distance);
                }
            }
        }
    }
}
