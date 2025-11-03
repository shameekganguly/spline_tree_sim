#ifndef COLLISION_CURSOR_H
#define COLLISION_CURSOR_H

#include "QuadraticSplineKinematic.h"
#include "SplineContact.h"
#include "TreeKinematic.h"

#include <Eigen/Core>

namespace spline_sim {

// struct to encapsulate response of getClosestBranchToPoint
struct CursorDistanceInfo {
  // branch name
  std::string branch_name;

  // contact point co-ordinates in spline frame
  SplinePointCartesian point_in_branch;

  // distance
  double distance;

  // ctor
  CursorDistanceInfo() : branch_name(""), distance(0.0) { /* Nothing to do */ }

  // ctor with parameters
  CursorDistanceInfo(std::string set_branch_name,
                     SplinePointCartesian set_point_in_branch,
                     double set_distance)
      : branch_name(set_branch_name), point_in_branch(set_point_in_branch),
        distance(set_distance) { /* Nothing to do */ }

  // print
  void print();
};

// Get closest branch to haptic cursor in a spline tree.
CursorDistanceInfo getClosestBranchToPoint(const TreeKinematic *tree,
                                           const Eigen::Vector3d &point);

} // namespace spline_sim

#endif
