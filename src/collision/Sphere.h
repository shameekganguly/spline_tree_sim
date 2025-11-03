#ifndef COLLISION_SPHERE_H
#define COLLISION_SPHERE_H

#include "QuadraticSplineKinematic.h"
#include "SplineContact.h"
#include "TreeKinematic.h"

#include <Eigen/Core>

namespace spline_sim {

// Sphere - TreeKinematic collision detection
void detectCollisionTreeSphere(std::vector<ContactInfo> &ret_list,
                               const TreeKinematic *tree,
                               const Eigen::Vector3d &sphere_center,
                               const double radius);

} // namespace spline_sim

#endif
