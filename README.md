# Ghostline Refitter

Rust library for fitting layered clothing in a shared character bind space.
There is no standalone application or game-resource decoder in this crate.

The character builder supplies body and garment geometry, native garment scores,
support paths and skeleton guides. Refitter computes geometric constraints and
displacements. The caller applies those changes to preview and export resources.

## Use the library

Add the crate to a sibling Rust project:

```toml
[dependencies]
ghostline-refitter = { path = "../refitter" }
```

The Rust import name is `ghostline_refitter`. Geometry uses `glam::Vec3`
in metres; `Triangle` describes three points in the common bind coordinate
system. Do not mix meshes in different coordinate frames.

Build the API reference from this repository:

```powershell
cargo doc --no-deps --open
```

| Module | Purpose |
| --- | --- |
| [surface](src/surface.rs) | Surface queries and collision geometry |
| [layers](src/layers.rs) | Clothing layer constraints |
| [envelope](src/envelope.rs) | Contacts measured in posed geometry |
| [armature](src/armature.rs) | Skeleton guidance |
| [intersection](src/intersection.rs) | Intersection measurements |
| [topology](src/topology.rs) | Mesh connectivity |

## Fitting model

The body remains a fixed collision boundary. The solver reserves space for lower
layers and compresses garments beneath higher-scored layers along their support
paths. A bounding-volume hierarchy accelerates surface queries.

The envelope module maps posed contacts back through calibrated skin transforms.
Callers combine representative poses, bound corrections in bind space and
inspect remaining contacts. Outward-normal checks distinguish surface contacts
from native closing caps.

Materials, UVs, indices, skin weights, garment flags, support buffers and physics
bindings remain the caller's responsibility. The library does not reproduce the
game's runtime garment-support simulation, retarget between skeletons or simulate
cloth. A successful solve does not certify animated clothing fit.

## Development

Run the independent geometry checks from this repository:

```powershell
cargo test --locked
```
