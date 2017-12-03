// TreeKinematic.h: Class for data encapsulation and kinematic computations
// for the whole tree

#ifndef TREE_KINEMATIC_H
#define TREE_KINEMATIC_H

#include "BranchKinematic.h"
#include <Eigen/Core>
#include <vector>
#include <map>
#include <string>

class TreeKinematic {
	// typedefs
	typedef std::map<std::string, BranchKinematic*> BranchList;
	typedef std::map<std::string, std::string> ParentNameMap;
	typedef std::map<std::string, Fruit*> FruitList;

	// public member functions
public:
	// ctor
	TreeKinematic(const std::string& name, const std::string& trunk_branch_name="trunk")
	: _name(name)
	{
		// create trunk branch
		_trunk = new BranchKinematic(trunk_branch_name);
		_branches[trunk_branch_name] = _trunk;
	}

	// dtor
	// NOTE: deleting the tree deletes all branches and fruits CURRENTLY
	// managed by this tree
	virtual ~TreeKinematic() {
		for (auto b_it: _branches) {
			auto branch = b_it.second;
			b_it.second = NULL;
			delete branch;
		}

		for (auto f_it: _fruits) {
			auto fruit = f_it.second;
			f_it.second = NULL;
			delete fruit;
		}

		delete _trunk;
	}

	/* ---- Member access/augment ---- */
	// get name
	const std::string& name() const {
		return _name;
	}

	// add branch
	BranchKinematic* branchIs(const std::string& name, const std::string& parent_name, double s);

	// access branch
	BranchKinematic* branch(const std::string& name) {
		auto it = _branches.find(name);
		if (it != _branches.end()) {
			return it->second;
		} else {
			return NULL;
		}
	}

	// const access branch
	const BranchKinematic* branch(const std::string& name) const {
		auto it = _branches.find(name);
		if (it != _branches.end()) {
			return it->second;
		} else {
			return NULL;
		}
	}

	// get non-const iterators
	BranchList::iterator branchesItrBegin() {
		return _branches.begin();
	}
	BranchList::iterator branchesItrEnd() {
		return _branches.end();
	}

	// get const iterators
	BranchList::const_iterator branchesItrBegin() const {
		return _branches.cbegin();
	}
	BranchList::const_iterator branchesItrEnd() const {
		return _branches.cend();
	}

	// remove branch. removes and returns all children branches as well
	// WARNING: This can cause a memory leak!! Because the returned pointers
	// point to memory no longer managed by this tree.
	BranchList branchRem(const std::string& name);

	// get parent branch name
	std::string branchParent(const std::string& name) const {
		auto it = _parent_map.find(name);
		if (it != _parent_map.end()) {
			return it->second;
		} else {
			return "";
		}
	}

	// update branch info in parent
	void branchInfoInParentUpdate(const std::string& branch_name, const ChildBranchInfo& new_info) {
		// get parent for branch
		auto it = _parent_map.find(branch_name);
		if (it != _parent_map.end()) {
			// update info
			auto parent = branch(it->second);
			parent->childBranchInfoIs(branch_name, new_info);
		}
	}

	// add fruit
	Fruit* fruitIs(const std::string& name, const std::string& branch_name, double s);

	// access fruit
	Fruit* fruit(const std::string& name) {
		auto it = _fruits.find(name);
		if (it != _fruits.end()) {
			return it->second;
		} else {
			return NULL;
		}
	}

	// const access branch
	const Fruit* fruit(const std::string& name) const {
		auto it = _fruits.find(name);
		if (it != _fruits.end()) {
			return it->second;
		} else {
			return NULL;
		}
	}

	// get non-const iterators
	FruitList::iterator fruitsItrBegin() {
		return _fruits.begin();
	}
	FruitList::iterator fruitsItrEnd() {
		return _fruits.end();
	}

	// get const iterators
	FruitList::const_iterator fruitsItrBegin() const {
		return _fruits.cbegin();
	}
	FruitList::const_iterator fruitsItrEnd() const {
		return _fruits.cend();
	}

	// get parent branch name for fruit
	std::string fruitParentBranch(const std::string& name) const {
		auto it = _fruit_branch_map.find(name);
		if (it != _fruit_branch_map.end()) {
			return it->second;
		} else {
			return "";
		}
	}

	// update fruit info in parent branch
	void fruitInfoInParentBranchUpdate(const std::string& fruit_name, const FruitInfo& new_info) {
		// get parent for branch
		auto it = _fruit_branch_map.find(fruit_name);
		if (it != _fruit_branch_map.end()) {
			// update info
			auto br = branch(it->second);
			br->fruitInfoIs(fruit_name, new_info);
		}
	}

	// remove fruit
	// WARNING: This can cause a memory leak!! Because the returned pointer
	// points to memory no longer managed by this tree.
	Fruit* fruitRem(const std::string& name);

	// get trunk name
	const std::string& trunk() const {
		return _trunk->name();
	}

	// // get trunk
	// BranchKinematic* trunk() {
	// 	return _trunk;
	// }

	// // get trunk const
	// const BranchKinematic* trunk() const {
	// 	return _trunk;
	// }

	// set transform
	void transformFromWorldIs(const Eigen::Affine3d& transform) {
		_transform = transform;
	}

	// get transform
	const Eigen::Affine3d& transformFromWorld() const {
		return _transform;
	}

	/* ---- Kinematic evaluations ---- */
	// get position in tree frame
	virtual void positionInTree(Eigen::Vector3d& ret_vec, const std::string& branch_name, double s) const;

	// get orientation in tree frame
	virtual void orientationInTree(Eigen::Vector3d& ret_vec, const std::string& branch_name, double s) const;

	// get position in world frame
	virtual void positionInWorld(Eigen::Vector3d& ret_vec, const std::string& branch_name, double s) const;

	// get orientation in world frame
	virtual void orientationInWorld(Eigen::Vector3d& ret_vec, const std::string& branch_name, double s) const;

	// get linear jacobian
	virtual void jacobianLinear(Eigen::MatrixXd& ret_mat, const std::string& branch_name, double s) const;

	// data members
protected:
	// name of tree
	std::string _name;

	// base frame pose in world frame
	Eigen::Affine3d _transform;

	// flat list of all branches
	BranchList _branches;

	// map of parents for branches. excludes root
	ParentNameMap _parent_map;

	// map of branches for fruits.
	ParentNameMap _fruit_branch_map;

	// flat list of all fruits
	FruitList _fruits;

	// trunk branch name
	BranchKinematic* _trunk;
};

#endif //TREE_KINEMATIC_H
