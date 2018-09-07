use bounding_volume::BoundingVolume;
use math::{Isometry, Point};
use na::Real;
use partitioning::{Visitor, VisitStatus};
use query::PointQuery;
use shape::CompositeShape;

/// Visitor for checking if a composite shape contains a specific point.
pub struct PointContainmentTest<'a, N: 'a + Real, S: 'a + CompositeShape<N>> {
    /// The composite shape on which the point containment test should be performed.
    pub shape: &'a S,
    /// The point to be tested.
    pub point: &'a Point<N>,
    /// A traversal will set this to `true` if the point is inside of `self.shape`.
    pub found: bool,
}

impl<'a, N: Real, BV: BoundingVolume<N> + PointQuery<N>, S: CompositeShape<N>> Visitor<usize, BV> for PointContainmentTest<'a, N, S> {
    #[inline]
    fn visit(&mut self, bv: &BV, b: Option<&usize>) -> VisitStatus {
        if bv.contains_point(&Isometry::identity(), self.point) {
            if let Some(b) = b {
                self.shape.map_part_at(*b, &mut |_, objm, obj| {
                    if obj.contains_point(objm, self.point) {
                        self.found = true;
                    }
                })
            }

            if self.found {
                VisitStatus::ExitEarly
            } else {
                VisitStatus::Continue
            }
        } else {
            VisitStatus::Stop
        }
    }
}