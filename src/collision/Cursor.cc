#include "Cursor.h"

#include <iostream>

#include <Eigen/Core>

#include "QuadraticSplineKinematic.h"
#include "SplineContact.h"
#include "TreeKinematic.h"

namespace spline_sim {

void CursorDistanceInfo::print() {
  std::cout << "-- cursor distance info --" << "\n";
  std::cout << "branch name " << branch_name << "\n";
  std::cout << "point in branch: s " << point_in_branch.s << " py "
            << point_in_branch.py << " pz " << point_in_branch.pz << "\n";
  std::cout << "distance " << distance << std::endl;
}

CursorDistanceInfo getClosestBranchToPoint(const TreeKinematic *tree,
                                           const Eigen::Vector3d &point) {
  CursorDistanceInfo ret_info;
  TreeKinematic::BranchList::const_iterator br_itr;
  Eigen::Affine3d branch_frame;
  SplinePointCartesian ret_point;
  double distance;
  Eigen::Vector3d normal;
  for (br_itr = tree->branchesItrBegin(); br_itr != tree->branchesItrEnd();
       ++br_itr) {
    // get position of point in local branch co-ordinates
    // TODO: this can be made more efficient by getting the branch
    // transforms only locally
    tree->transformInWorld(branch_frame, br_itr->first, 0.0);
    br_itr->second->spline()->closestPointToPoint(
        ret_point, distance, normal, branch_frame.inverse() * point);
    // check for minimum distance
    if (ret_info.branch_name.empty() || distance < ret_info.distance) {
      ret_info.branch_name = br_itr->first;
      ret_info.point_in_branch = ret_point;
      ret_info.distance = distance;
    }
  }
  return ret_info;
}

} // namespace spline_sim