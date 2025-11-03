// TreeVisual.h: Implementation class to manage the visualization of the
// branches and fruits on a given tree

#ifndef TREE_VISUAL_H
#define TREE_VISUAL_H

#include <string>

#include <chai3d.h>

#include "TreeKinematic.h"

namespace spline_sim {

// NOTE: we do not maintain a tree list on chai. we only maintain one in
// here.
// Also, this class resets its local transform with respect to the world
// always to identity because we maintain a transform on the kinematic
// element.
class TreeVisual : public chai3d::cGenericObject {
public:
  // ctor
  TreeVisual(TreeKinematic *kinematic)
      : _kinematic(kinematic) { /* Nothing to do */ }

  // update graphics
  // updates the position and orientation on all elements
  void updateGraphics();

  // set material for all branches. convenience function
  void branchMaterialIs(chai3d::cMaterialPtr material);

  // set material for all fruits. convenience function
  void fruitMaterialIs(chai3d::cMaterialPtr material);

public:
  // kinematic tree
  TreeKinematic *_kinematic;
};

} // namespace spline_sim

#endif // TREE_VISUAL_H