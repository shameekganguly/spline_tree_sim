#ifndef TREE_DYNAMICS_H
#define TREE_DYNAMICS_H

#include <optional>

#include "LCPSolver.h"
#include "TreeKinematic.h"
#include "collision/SplineContact.h"

namespace spline_sim {

class TreeDynamics {
public:
  typedef std::pair<std::string, Eigen::Vector3d> FruitHapticForce;

  TreeDynamics(TreeKinematic *tree_kinematic)
      : tree_kinematic_(tree_kinematic) {}

  void Step(double dt_sec, const std::vector<ContactInfo> &contact_list,
            const Eigen::VectorXd &gamma_external);

  void Step(double dt_sec, const std::vector<ContactInfo> &contact_list);

  void Step(double dt_sec, const std::vector<ContactInfo> &contact_list,
            const FruitHapticForce &fruit_haptic_force);

protected:
  TreeKinematic *tree_kinematic_;
  LCPSolver contact_solver;
};

} // namespace spline_sim

#endif
