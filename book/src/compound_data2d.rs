extern crate nalgebra;
extern crate "ncollide2df32" as ncollide;

use std::sync::Arc;
use nalgebra::na::{Iso2, Vec2};
use nalgebra::na;
use ncollide::geom::{Geom, Plane, Cuboid, Compound, CompoundData};
use ncollide::volumetric::Volumetric;

fn main() {
    // delta transformation matrices.
    let delta1 = Iso2::new(Vec2::new(0.0f32, -1.5), na::zero());
    let delta2 = Iso2::new(Vec2::new(-1.5f32, 0.0), na::zero());
    let delta3 = Iso2::new(Vec2::new(1.5f32,  0.0), na::zero());

    // 1) initialize the CompoundData.
    let mut compound_data = CompoundData::new();

    /*
     * push_geom
     */
    // The mass, center of mass and angular inertia tensor are automatically
    // computed.
    compound_data.push_geom(delta1, Cuboid::new(Vec2::new(1.5f32, 0.75)), 1.0);

    /*
     * push_geom_with_mass_properties
     */
    // area                   = 1.0
    // mass                   = 10.0
    // center of mass         = the origin (na::zero())
    // angular inertia tensor = identity matrix (na::one())
    compound_data.push_geom_with_mass_properties(
        delta2,
        Plane::new(Vec2::new(1f32, 0.0)),
        (na::one(), 10.0f32, na::one(), na::one()));

    /*
     * push_shared_geom_with_mass_properties
     */
    // The geometry we want to share.
    let cuboid = Cuboid::new(Vec2::new(0.75f32, 1.5));
    // Make ncollide compute the mass properties of the cuboid.
    let (c_mass, c_com, c_tensor) = cuboid.mass_properties(&1.0); // density = 1.0
    let c_area                    = cuboid.surface();
    // Build the shared geometry.
    let shared_cuboid = Arc::new(box cuboid as Box<Geom + Send + Sync>);
    // Add the geometry to the compound data.
    compound_data.push_shared_geom_with_mass_properties(
        delta3,
        shared_cuboid.clone(),
        (c_area, c_mass, c_com, c_tensor));
    // `shared_cuboid` can still be used thereafter…

    // 2) create the compound geometry.
    let compound = Compound::new(compound_data);

    assert!(compound.geoms().len() == 3);
}
