#ifndef COLLISION_PLANE_H
#define COLLISION_PLANE_H

#include "QuadraticSplineKinematic.h"
#include "SplineContact.h"
#include "TreeKinematic.h"

#include <Eigen/Core>

namespace spline_sim {

// Plane - TreeKinematic collision detection
void detectCollisionTreePlane(std::vector<ContactInfo> &ret_list,
                               const TreeKinematic *tree,
                               const Eigen::Vector3d &plane_point,
                               const Eigen::Vector3d &plane_normal);

} // namespace spline_sim

#endif
