#include "Sphere.h"

#include <limits>

#include "QuadraticSplineKinematic.h"
#include "SplineContact.h"
#include "TreeKinematic.h"

#include <Eigen/Core>

using namespace Eigen;

namespace spline_sim {

namespace {
constexpr double kDistanceTol = 1e-3;

ContactInfo makeContactInfo(const Eigen::Vector3d &plane_normal_in_spline,
                            const QuadraticSplineKinematic &spline,
                            double spline_pt_s,
                            double plane_pt_normal_distance) {
  Matrix3d orientation_s;
  spline.splineOrientation(orientation_s, spline_pt_s);
  Vector3d normal_local = orientation_s.transpose() * plane_normal_in_spline;

  // Note: normal_local[0] will be zero if spline_pt_s lies in [0, _length]
  SplinePointCartesian point_in_spline;
  point_in_spline.s = spline_pt_s - spline._radius * normal_local[0];
  point_in_spline.py = -spline._radius * normal_local[1];
  point_in_spline.pz = -spline._radius * normal_local[2];

  ContactInfo info;
  info.point_in_branch = point_in_spline;
  info.normal = -plane_normal_in_spline;
  double distance = plane_pt_normal_distance - spline._radius;
  info.penetration_depth = (distance < 0) ? -distance : 0;
  return info;
}

bool detectCollisionSplinePlane3points(
    std::vector<ContactInfo> &ret_list, const QuadraticSplineKinematic &spline,
    const Eigen::Vector3d &spline_mid_pt, const Eigen::Vector3d &spline_end_pt,
    const Eigen::Vector3d &plane_point_in_spline,
    const Eigen::Vector3d &plane_normal_in_spline) {
  // The origin, mid-pt and end-pt should be equidistant from the plane
  double origin_dist = -plane_point_in_spline.dot(plane_normal_in_spline);
  double mid_pt_dist =
      (spline_mid_pt - plane_point_in_spline).dot(plane_normal_in_spline);
  double end_pt_dist =
      (spline_end_pt - plane_point_in_spline).dot(plane_normal_in_spline);
  // std::cout << "detectCollisionSplinePlane3points: origin_dist " <<
  // origin_dist
  // << " mid_pt_dist " << mid_pt_dist << " end_pt_dist " << end_pt_dist
  // << "\n";
  if ((origin_dist < spline._radius + kDistanceTol) &&
      (mid_pt_dist < spline._radius + kDistanceTol) &&
      (end_pt_dist < spline._radius + kDistanceTol)) {
    ret_list.push_back(makeContactInfo(plane_normal_in_spline, spline,
                                       /*spline_pt_s=*/0, origin_dist));
    ret_list.push_back(makeContactInfo(plane_normal_in_spline, spline,
                                       /*spline_pt_s=*/0.5 * spline._length,
                                       mid_pt_dist));
    ret_list.push_back(makeContactInfo(plane_normal_in_spline, spline,
                                       /*spline_pt_s=*/spline._length,
                                       end_pt_dist));
    return true;
  }
  return false;
}

// Should be called only if detectCollisionSplinePlane3points returns false
bool detectCollisionSplinePlane2points(
    std::vector<ContactInfo> &ret_list, const QuadraticSplineKinematic &spline,
    const Eigen::Vector3d &spline_end_pt,
    const Eigen::Vector3d &plane_point_in_spline,
    const Eigen::Vector3d &plane_normal_in_spline) {
  // The origin, end-pt should be equidistant from the plane
  double origin_dist = -plane_point_in_spline.dot(plane_normal_in_spline);
  double end_pt_dist =
      (spline_end_pt - plane_point_in_spline).dot(plane_normal_in_spline);
  if ((origin_dist < spline._radius + kDistanceTol) &&
      (end_pt_dist < spline._radius + kDistanceTol)) {
    ret_list.push_back(makeContactInfo(plane_normal_in_spline, spline,
                                       /*spline_pt_s=*/0, origin_dist));
    ret_list.push_back(makeContactInfo(plane_normal_in_spline, spline,
                                       /*spline_pt_s=*/spline._length,
                                       end_pt_dist));
    return true;
  }
  return false;
}

// Should be called only if detectCollisionSplinePlane2points returns false
bool detectCollisionSplinePlane1point(
    std::vector<ContactInfo> &ret_list, const QuadraticSplineKinematic &spline,
    const Eigen::Vector3d &plane_point_in_spline,
    const Eigen::Vector3d &plane_normal_in_spline) {
  std::vector<double> test_s;
  test_s.push_back(0);
  test_s.push_back(spline._length);

  double gamma = spline.gam();
  // If gamma is zero, the spline is in its unmodified state. So either
  // `detectCollisionSplinePlane3points` returned true or one of the end points
  // is closest to the plane.
  if (fabs(gamma) > 1e-5) {

    // Solve spline_x_dir(s).dot(plane_normal_in_spline) = 0 for s
    double tb = tan(spline._beta);
    double sa = sin(spline._alpha);
    double t_den = sqrt(pow(tb, 2) + pow(sa, 2));
    double cos_t = tb / t_den;
    double sin_t = -sa / t_den;
    double num = -plane_normal_in_spline(0);
    double den =
        plane_normal_in_spline.segment<2>(1).dot(Vector2d(cos_t, sin_t));

    // If num and den are both zero, that indicates the test plane is
    // parallel to the spline plane. This should be covered in
    // `detectCollisionSplinePlane3points`, so we can return early.
    if (fabs(num) < 1e-5 && fabs(den) < 1e-5) {
      return false;
    }
    double phi = atan(num / den);
    // std::cout << "phi: " << phi << " gamma: " << gamma << "\n";
    // Both phi and phi + M_PI can be valid solutions.
    test_s.push_back(phi / gamma * spline._length);
    test_s.push_back((phi + M_PI) / gamma * spline._length);
  }

  double min_dist = std::numeric_limits<double>::max();
  double min_s;
  for (uint i = 0; i < test_s.size(); i++) {
    if (test_s[i] < 0 || test_s[i] > spline._length) {
      // std::cout << "Skip test_s: " << test_s[i] << "\n";
      continue;
    }
    Vector3d spline_pt;
    spline.splineLocation(spline_pt, test_s[i]);
    double dist =
        (spline_pt - plane_point_in_spline).dot(plane_normal_in_spline);
    // std::cout << "test_s : " << test_s[i] << " dist: " << dist << "\n";
    if (dist < min_dist) {
      min_dist = dist;
      min_s = test_s[i];
    }
  }
  if (min_dist < spline._radius + kDistanceTol) {
    ret_list.push_back(
        makeContactInfo(plane_normal_in_spline, spline, min_s, min_dist));
    return true;
  }
  return false;
}

} // namespace

void detectCollisionTreePlane(std::vector<ContactInfo> &ret_list,
                              const TreeKinematic *tree,
                              const Eigen::Vector3d &plane_point,
                              const Eigen::Vector3d &plane_normal) {
  ret_list.clear();
  for (auto br_itr = tree->branchesItrBegin(); br_itr != tree->branchesItrEnd();
       ++br_itr) {
    // get position of sphere_center in local branch co-ordinates
    // TODO: this can be made more efficient by getting the branch
    // transforms only locally
    Eigen::Affine3d branch_frame;
    tree->transformInWorld(branch_frame, br_itr->first, 0.0);
    Eigen::Vector3d plane_normal_in_spline =
        branch_frame.rotation().transpose() * plane_normal;
    Eigen::Vector3d plane_point_in_spline =
        branch_frame.inverse() * plane_point;

    const auto *spline_ptr = br_itr->second->spline();

    // Broadphase: check if mid-point is > (R + L/2) distance away from plane
    Vector3d mid_point;
    spline_ptr->splineLocation(mid_point, 0.5 * spline_ptr->_length);
    if (fabs((mid_point - plane_point_in_spline).dot(plane_normal_in_spline)) >
        (spline_ptr->_radius + spline_ptr->_length / 2)) {
      continue;
    }

    Vector3d end_point;
    spline_ptr->splineLocation(end_point, spline_ptr->_length);

    std::vector<ContactInfo> branch_ret_list;
    if (detectCollisionSplinePlane3points(
            branch_ret_list, *spline_ptr, mid_point, end_point,
            plane_point_in_spline, plane_normal_in_spline)) {
      // std::cout << "Here1\n";
    } else if (detectCollisionSplinePlane2points(
                   branch_ret_list, *spline_ptr, end_point,
                   plane_point_in_spline, plane_normal_in_spline)) {
      // std::cout << "Here2\n";
    } else if (detectCollisionSplinePlane1point(branch_ret_list, *spline_ptr,
                                                plane_point_in_spline,
                                                plane_normal_in_spline)) {
      // std::cout << "Here3\n";
    }

    for (uint i = 0; i < branch_ret_list.size(); i++) {
      branch_ret_list[i].normal =
          branch_frame.rotation() * branch_ret_list[i].normal;
      branch_ret_list[i].branch_name = br_itr->first;
      ret_list.push_back(branch_ret_list[i]);
    }
  }
}

} // namespace spline_sim
