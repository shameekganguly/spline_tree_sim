#ifndef SIMULATION_WORLD_H
#define SIMULATION_WORLD_H

#include <Eigen/Core>

#include "TreeDynamics.h"
#include "TreeKinematic.h"
#include "collision/SplineContact.h"

namespace spline_sim {
struct PlaneInfo {
  std::string name;
  Eigen::Vector3d point;
  Eigen::Vector3d normal;
};

class World {
public:
  typedef std::unordered_map<std::string, std::unique_ptr<TreeKinematic>>
      TreeKinematicMap;

  void Collide();
  void Step(double dt_sec);

  // Takes ownership
  void AddTree(std::unique_ptr<TreeKinematic> tree);

  void AddPlane(PlaneInfo plane);

  TreeKinematicMap::iterator TreesItrBegin() {
    return trees_kinematic_.begin();
  }
  TreeKinematicMap::iterator TreesItrEnd() { return trees_kinematic_.end(); }

  TreeKinematicMap::const_iterator TreesItrBegin() const {
    return trees_kinematic_.cbegin();
  }
  TreeKinematicMap::const_iterator TreesItrEnd() const {
    return trees_kinematic_.cend();
  }

protected:
  TreeKinematicMap trees_kinematic_;
  std::unordered_map<std::string, std::unique_ptr<TreeDynamics>> trees_dynamic_;
  std::unordered_map<std::string, PlaneInfo> planes_;

  // Updated by calling `Collide`. Passed to `tree_dynamic_->Step()`.
  std::unordered_map<std::string, std::vector<ContactInfo>> contact_lists_;
};
} // namespace spline_sim

#endif