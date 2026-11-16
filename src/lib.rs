//! Reusable fitting constraints for native garments sharing a bind-space skeleton.
//! The caller retains materials, UVs, indices, skinning, support and physics data.

pub mod armature;
pub mod envelope;
pub mod intersection;
pub mod layers;
pub mod surface;
pub mod topology;

use glam::Vec3;

/// A triangle in the common garment bind coordinate system, in metres.
#[derive(Clone, Copy, Debug)]
pub struct Triangle(pub [Vec3; 3]);

/// Measured geometry changes from a fitting pass.
#[derive(Default, Debug)]
pub struct FitReport {
    pub changed_vertices: usize,
    pub maximum_displacement: f32,
}
