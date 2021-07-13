use na::{RealField, Unit};

use crate::interpolation::RigidMotion;
use crate::math::{Isometry, Point, Vector};
use crate::query::{self, ClosestPoints, TOIStatus, TOI};
use crate::shape::SupportMap;

/// Time of impacts between two support-mapped shapes under a rigid motion.
pub fn nonlinear_time_of_impact_support_map_support_map<N, G1: ?Sized, G2: ?Sized>(
    motion1: &(impl RigidMotion<N> + ?Sized),
    g1: &G1,
    motion2: &(impl RigidMotion<N> + ?Sized),
    g2: &G2,
    max_toi: N,
    target_distance: N,
) -> Option<TOI<N>>
where
    N: RealField,
    G1: SupportMap<N>,
    G2: SupportMap<N>,
{
    nonlinear_time_of_impact_support_map_support_map_with_closest_points_function(
        motion1,
        g1,
        motion2,
        g2,
        max_toi,
        target_distance,
        query::closest_points_support_map_support_map,
    )
}

/// Time of impacts between two support-mapped shapes under a rigid motion.
///
/// You probably want to use `query::nonlinear_time_of_impact_support_map_support_map` instead of this one.
/// The distance function between the two shapes must be given.
pub fn nonlinear_time_of_impact_support_map_support_map_with_closest_points_function<
    N,
    G1: ?Sized,
    G2: ?Sized,
>(
    motion1: &(impl RigidMotion<N> + ?Sized),
    g1: &G1,
    motion2: &(impl RigidMotion<N> + ?Sized),
    g2: &G2,
    max_toi: N,
    target_distance: N,
    closest_points: impl Fn(&Isometry<N>, &G1, &Isometry<N>, &G2, N) -> ClosestPoints<N>,
) -> Option<TOI<N>>
where
    N: RealField,
    G1: SupportMap<N>,
    G2: SupportMap<N>,
{
    let _0_5 = na::convert(0.5);
    let mut min_t = N::zero();
    let mut prev_min_t = min_t;
    let abs_tol: N = query::algorithms::gjk::eps_tol();
    let rel_tol = abs_tol.sqrt();
    let mut result = TOI {
        toi: N::zero(),
        normal1: Vector::x_axis(),
        normal2: Vector::x_axis(),
        witness1: Point::origin(),
        witness2: Point::origin(),
        status: TOIStatus::Penetrating,
    };

    loop {
        let pos1 = motion1.position_at_time(result.toi);
        let pos2 = motion2.position_at_time(result.toi);

        // FIXME: use the _with_params version of the closest points query.
        match closest_points(&pos1, g1, &pos2, g2, N::max_value()) {
            ClosestPoints::Intersecting => {
                if result.toi == N::zero() {
                    result.status = TOIStatus::Penetrating
                } else {
                    result.status = TOIStatus::Failed;
                }
                break;
            }
            ClosestPoints::WithinMargin(p1, p2) => {
                // FIXME: do the "inverse_transform_point" only when we are about to return.
                // the result.
                result.witness1 = pos1.inverse_transform_point(&p1);
                result.witness2 = pos2.inverse_transform_point(&p2);

                if let Some((dir, mut dist)) = Unit::try_new_and_get(p2 - p1, N::default_epsilon())
                {
                    // FIXME: do the "inverse transform unit vector" only when we are about to return.
                    result.normal1 = pos1.inverse_transform_unit_vector(&dir);
                    result.normal2 = pos2.inverse_transform_unit_vector(&-dir);

                    let mut niter = 0;
                    min_t = result.toi;
                    let mut max_t = max_toi;
                    let min_target_distance = (target_distance - rel_tol).max(N::zero());
                    let max_target_distance = target_distance + rel_tol;

                    loop {
                        // FIXME: use the secant method too for finding the next iterate.
                        if dist < min_target_distance {
                            // Too close or penetration, go back in time.
                            max_t = result.toi;
                            result.toi = (min_t + result.toi) * _0_5;
                        } else if dist > max_target_distance {
                            // Too far apart, go forward in time.
                            min_t = result.toi;
                            result.toi = (result.toi + max_t) * _0_5;
                        } else {
                            // Reached tolerance, break.
                            break;
                        }

                        if max_t - min_t < max_t.max(na::convert(1.)) * abs_tol {
                            result.toi = min_t;
                            break;
                        }

                        let pos1 = motion1.position_at_time(result.toi);
                        let pos2 = motion2.position_at_time(result.toi);
                        let pt1 = g1.support_point_toward(&pos1, &dir);
                        let pt2 = g2.support_point_toward(&pos2, &-dir);

                        dist = (pt2 - pt1).dot(&dir);

                        niter += 1;
                    }

                    if min_t - prev_min_t < abs_tol {
                        if max_t == max_toi {
                            let pos1 = motion1.position_at_time(max_t);
                            let pos2 = motion2.position_at_time(max_t);
                            // Check the configuration at max_t to see if the object are not disjoint.
                            // NOTE: could we do this earlier, before the above loop?
                            // It feels like this could prevent catching some corner-cases like
                            // if one object is rotated by almost 180 degrees while the other is immobile.
                            let pt1 = g1.support_point_toward(&pos1, &dir);
                            let pt2 = g2.support_point_toward(&pos2, &-dir);
                            if (pt2 - pt1).dot(&dir) > target_distance {
                                // We found an axis that separate both objects at the end configuration.
                                return None;
                            }
                        }

                        result.status = TOIStatus::Converged;
                        break;
                    }

                    prev_min_t = min_t;

                    if niter == 0 {
                        result.status = TOIStatus::Converged;
                        break;
                    }
                } else {
                    result.status = TOIStatus::Failed;
                    break;
                }
            }
            ClosestPoints::Disjoint => unreachable!(),
        }
    }

    Some(result)
}
