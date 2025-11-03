// A utility file for spline contact implementations
#ifndef SPLINE_CONTACT_H
#define SPLINE_CONTACT_H

#include "QuadraticSplineKinematic.h"
#include "TreeKinematic.h"

#include <Eigen/Core>
#include <string>

namespace spline_sim {

// contact between haptic sphere cursor and a branch on the tree
struct ContactInfo {
  // branch name
  std::string branch_name;

  // contact point co-ordinates in spline frame
  SplinePointCartesian point_in_branch;

  // outward normal at contact point in the global frame
  Eigen::Vector3d normal;

  // penetration depth
  double penetration_depth;

  // ctor
  ContactInfo()
      : branch_name(""), penetration_depth(0.0) { /* Nothing to do */ }

  // ctor with parameters
  ContactInfo(std::string set_branch_name,
              SplinePointCartesian set_point_in_branch,
              double set_penetration_depth)
      : branch_name(set_branch_name), point_in_branch(set_point_in_branch),
        penetration_depth(set_penetration_depth) { /* Nothing to do */ }

  // print
  void print();
};

} // namespace spline_sim

#endif // SPLINE_CONTACT_H
