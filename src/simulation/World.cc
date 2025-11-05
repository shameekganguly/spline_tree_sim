#include "World.h"
#include "collision/Plane.h"

namespace spline_sim {

void World::Collide() {
  for (auto tree_itr = trees_kinematic_.begin();
       tree_itr != trees_kinematic_.end(); tree_itr++) {
    const std::string &tree_name = tree_itr->first;
    auto tree_ptr = tree_itr->second.get();
    if (contact_lists_.find(tree_name) == contact_lists_.end()) {
      contact_lists_[tree_name] = {};
    } else {
      contact_lists_[tree_name].clear();
    }
    for (auto plane_itr = planes_.begin(); plane_itr != planes_.end();
         plane_itr++) {
      std::vector<ContactInfo> contacts;
      detectCollisionTreePlane(contacts, tree_ptr, plane_itr->second.point,
                               plane_itr->second.normal);
      for (auto contact : contacts) {
        contact_lists_[tree_name].push_back(contact);
      }
    }
  }
}

void World::Step(double dt_sec) {
  for (auto tree_itr = trees_dynamic_.begin(); tree_itr != trees_dynamic_.end();
       tree_itr++) {
    tree_itr->second->Step(dt_sec, contact_lists_[tree_itr->first]);
  }
}

void World::AddTree(std::unique_ptr<TreeKinematic> tree) {
  trees_dynamic_.emplace(tree->name(),
                         std::make_unique<TreeDynamics>(tree.get()));
  trees_kinematic_.emplace(tree->name(), std::move(tree));
}

void World::AddPlane(PlaneInfo plane) { planes_.emplace(plane.name, plane); }

} // namespace spline_sim
