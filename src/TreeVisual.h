// TreeVisual.h: Implementation class to manage the visualization of the 
// branches and fruits on a given tree

#ifndef TREE_VISUAL_H
#define TREE_VISUAL_H

#include "TreeKinematic.h"
#include <world/CGenericObject.h>

#include <Eigen/Core>
#include <deque>
#include <string>

#include <iostream>
using namespace std;
using namespace chai3d;
using namespace Eigen;

// NOTE: we do not maintain a tree list on chai. we only maintain one in
// here.
// Also, this class resets its local transform with respect to the world
// always to identity because we maintain a transform on the kinematic
// element.
class TreeVisual: public cGenericObject {
public:
	// ctor
	TreeVisual(TreeKinematic* kinematic)
	: _kinematic(kinematic)
	{ /* Nothing to do */ }

	// update graphics
	// updates the position and orientation on all elements
	void updateGraphics(cGenericObject* chai_parent) {
		// reset local position and orientation if it has been changed
		setLocalPos(Vector3d::Zero());
		setLocalRot(Matrix3d::Identity());

		// cache variables
		std::deque<std::string> list_to_process;
		BranchKinematic* br;
		Fruit* fr;
		BranchKinematic* parent_br;
		QuadraticSplineVisual* br_visual;
		Vector3d position, fruit_position;
		Matrix3d rotation, fruit_rotation;

		// process trunk
		list_to_process.push_back(_kinematic->trunk());

		// process rest
		while (!list_to_process.empty()) {
			// get branch from front of queue
			string br_name = list_to_process.front();
			br = _kinematic->branch(br_name);
			list_to_process.pop_front();

			// compute the branch transform
			br_visual = br->splineVisual();
			if (br_name.compare(_kinematic->trunk()) == 0) {
				// process trunk separately
				position = _kinematic->transformFromWorld().translation();
				rotation = _kinematic->transformFromWorld().rotation();
			} else {
				parent_br = _kinematic->branch(_kinematic->branchParent(br_name));
				ChildBranchInfo child_info = parent_br->childBranchInfo(br_name);
				parent_br->spline()->splineOrientation(rotation, child_info.s);
				rotation = parent_br->splineVisual()->getLocalRot().eigen() * rotation;
				parent_br->spline()->splineLocation(position, child_info.s);
				position = rotation*position;
				position += parent_br->splineVisual()->getLocalPos().eigen();
				rotation = rotation * child_info.rotation;
			}
			br_visual->setLocalPos(position);
			br_visual->setLocalRot(rotation);

			// add branch to chai parent object if it does not exist already
			bool add_to_parent = false;
			auto br_visual_parent = br_visual->getParent();
			if (br_visual_parent) {
				if (br_visual_parent->m_name.compare(chai_parent->m_name) != 0) {
					br_visual_parent->removeChild(br_visual);
					add_to_parent = true;
				}
			} else { add_to_parent = true; }
			if (add_to_parent) {
				chai_parent->addChild(br_visual);
			}

			// compute transforms for fruits
			for (auto fruit_info_itr: br->_fruits) {
				FruitInfo fruit_info = fruit_info_itr.second;
				fr = _kinematic->fruit(fruit_info.name);
				auto fruit_visual = fr->graphic();
				br->spline()->splineOrientation(fruit_rotation, fruit_info.s);
				fruit_rotation = rotation * fruit_rotation;
				br->spline()->splineLocation(fruit_position, fruit_info.s);
				fruit_position = position + rotation * fruit_position;
				fruit_visual->setLocalRot(fruit_rotation);
				fruit_visual->setLocalPos(fruit_position);

				// add fruit to chai parent if it does not exist already
				bool add_fruit_to_parent = false;
				auto fruit_visual_parent = fruit_visual->getParent();
				if (fruit_visual_parent) {
					if (fruit_visual_parent->m_name.compare(chai_parent->m_name) != 0) {
						fruit_visual_parent->removeChild(fruit_visual);
						add_fruit_to_parent = true;
					}
				} else { add_fruit_to_parent = true; }
				if (add_fruit_to_parent) {
					chai_parent->addChild(fruit_visual);
				}
			}

			// append children to list
			for (auto it_ch: br->_children) {
				auto child_name = it_ch.first;
				list_to_process.push_back(child_name);
			}
		}
	}

	// set material for all branches. convenience function
	void branchMaterialIs(cMaterialPtr material) {
		// iterate over all branches
		for (auto itr=_kinematic->branchesItrBegin(); itr!=_kinematic->branchesItrEnd(); ++itr) {
			auto visual = itr->second->splineVisual();
			visual->m_material = material;
		}
	}

	// set material for all fruits. convenience function
	void fruitMaterialIs(cMaterialPtr material) {
		// iterate over all fruits
		for (auto itr=_kinematic->fruitsItrBegin(); itr!=_kinematic->fruitsItrEnd(); ++itr) {
			auto visual = itr->second->graphic();
			visual->m_material = material;
		}
	}

public:
	// kinematic tree
	TreeKinematic* _kinematic;
};

#endif //TREE_VISUAL_H