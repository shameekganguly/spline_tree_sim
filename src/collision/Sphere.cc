#include "Sphere.h"

#include "QuadraticSplineKinematic.h"
#include "SplineContact.h"
#include "TreeKinematic.h"

#include <Eigen/Core>

namespace spline_sim {

// Sphere - TreeKinematic collision detection
void detectCollisionTreeSphere(std::vector<ContactInfo> &ret_list,
                               const TreeKinematic *tree,
                               const Eigen::Vector3d &sphere_center,
                               const double radius) {
  ret_list.clear();
  TreeKinematic::BranchList::const_iterator br_itr;
  Eigen::Affine3d branch_frame;
  SplinePointCartesian ret_point;
  double distance;
  Eigen::Vector3d normal;
  for (br_itr = tree->branchesItrBegin(); br_itr != tree->branchesItrEnd();
       ++br_itr) {
    // get position of sphere_center in local branch co-ordinates
    // TODO: this can be made more efficient by getting the branch
    // transforms only locally
    tree->transformInWorld(branch_frame, br_itr->first, 0.0);
    br_itr->second->spline()->closestPointToSphere(
        ret_point, distance, normal, branch_frame.inverse() * sphere_center,
        radius);
    // check for collision
    if (distance < 1e-3) {
      ContactInfo info;
      info.branch_name = br_itr->first;
      info.point_in_branch = ret_point;
      info.normal = branch_frame.rotation() * normal;
      info.penetration_depth = (distance < 0) ? -distance : 0.0;
      ret_list.push_back(info);
    }
  }
}

} // namespace spline_sim
