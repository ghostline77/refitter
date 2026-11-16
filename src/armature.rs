//! Rigid rebasing of mechanical child bones onto named anatomical anchors.
use glam::Mat4;
use std::collections::BTreeMap;

/// A malformed hierarchy, anchor, or transform supplied to an armature fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArmatureError {
    /// Counts, names, parents or matrices fail validation.
    InvalidInput,
    /// Parent links contain a cycle.
    CyclicHierarchy,
}
impl std::fmt::Display for ArmatureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => "invalid armature fit input",
            Self::CyclicHierarchy => "cyclic armature hierarchy",
        })
    }
}
impl std::error::Error for ArmatureError {}

/// Replaces explicitly named model-space anchors and preserves every other
/// bone's local transform. Mechanical descendants therefore move rigidly with
/// their anatomical mounting joint rather than stretching to fit a body mesh.
/// All matrices must use the same coordinate system. Bone order need not be topological.
///
/// # Errors
/// Returns an error for invalid matrices, duplicate names, unmatched anchors,
/// out-of-range parents or a cyclic hierarchy.
pub fn rebase(
    names: &[String],
    parents: &[Option<usize>],
    models: &[Mat4],
    anchors: &BTreeMap<String, Mat4>,
) -> Result<Vec<Mat4>, ArmatureError> {
    let valid = |m: &Mat4| {
        m.is_finite()
            && m.inverse().is_finite()
            && (m.w_axis.w - 1.).abs() < 1e-5
            && m.x_axis.w.abs() < 1e-5
            && m.y_axis.w.abs() < 1e-5
            && m.z_axis.w.abs() < 1e-5
    };
    if names.is_empty()
        || names.len() != parents.len()
        || names.len() != models.len()
        || models.iter().chain(anchors.values()).any(|m| !valid(m))
        || parents.iter().flatten().any(|&p| p >= names.len())
        || names
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
            .len()
            != names.len()
        || anchors.is_empty()
        || anchors.keys().any(|name| !names.contains(name))
    {
        return Err(ArmatureError::InvalidInput);
    }
    let mut states = vec![0; names.len()];
    let mut result = vec![Mat4::IDENTITY; names.len()];
    for i in 0..names.len() {
        let mut chain = Vec::new();
        let mut current = Some(i);
        while let Some(j) = current {
            if states[j] == 2 {
                break;
            }
            if states[j] == 1 {
                return Err(ArmatureError::CyclicHierarchy);
            }
            states[j] = 1;
            chain.push(j);
            current = parents[j];
        }
        for j in chain.into_iter().rev() {
            result[j] = anchors.get(&names[j]).copied().unwrap_or_else(|| {
                parents[j].map_or(models[j], |p| result[p] * models[p].inverse() * models[j])
            });
            states[j] = 2;
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::{Quat, Vec3};
    #[test]
    fn a_rotated_mount_preserves_its_mechanism_and_the_other_arm() {
        let names = ["Root", "LeftForeArm", "hinge", "blade", "RightForeArm"].map(str::to_owned);
        let parents = [None, Some(0), Some(1), Some(2), Some(0)];
        let source = [
            Vec3::ZERO,
            Vec3::X,
            Vec3::new(1., 0.1, 0.),
            Vec3::new(1., 0.4, 0.),
            -Vec3::X,
        ]
        .map(Mat4::from_translation);
        let anchor =
            Mat4::from_rotation_translation(Quat::from_rotation_z(0.4), Vec3::new(1.3, 0.2, 0.));
        let fitted = rebase(
            &names,
            &parents,
            &source,
            &BTreeMap::from([("LeftForeArm".into(), anchor)]),
        )
        .unwrap();
        assert!(fitted[1].abs_diff_eq(anchor, 1e-6));
        assert!(
            (fitted[2].inverse() * fitted[3]).abs_diff_eq(source[2].inverse() * source[3], 1e-6)
        );
        assert_eq!(fitted[4], source[4]);
        assert!(
            rebase(
                &names,
                &[Some(3), Some(0), Some(1), Some(2), Some(0)],
                &source,
                &BTreeMap::from([("LeftForeArm".into(), anchor)])
            )
            .is_err()
        );
    }
}
