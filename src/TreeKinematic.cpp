#include "TreeKinematic.h"

#include <stdexcept>

#include <iostream>
using namespace std;
using namespace Eigen;

// add branch
BranchKinematic* TreeKinematic::branchIs(const std::string& name, const std::string& parent_name, double s) {
	// check if branch exists already.
	// if so, check parent
	// if parent is correct, simply return the branch
	// else error out
	auto it_p = _parent_map.find(name);
	if (it_p != _parent_map.end()) {
		if (parent_name.compare(it_p->second) == 0) {
			return _branches[name];
		} else {
			throw(std::runtime_error("Branch exists with different parent."));
		}
	}

	// check if parent exists.
	auto it_b = _branches.find(parent_name);
	if (it_b == _branches.end()) {
		throw(std::runtime_error("Parent branch does not exist."));
	}

	// check name conflict with trunk
	if (name.compare(trunk()) == 0) {
		throw(std::runtime_error("Cannot have same name as trunk."));
	}

	// else add branch
	BranchKinematic* new_branch = new BranchKinematic(name);
	// - create child branch info 
	ChildBranchInfo info;
	info.name = name;
	info.s = s;
	// - add to parent branch
	auto parent = branch(parent_name);
	parent->_children[name] = info;
	// - add branch to self
	_branches[name] = new_branch;
	// - update parent map
	_parent_map[name] = parent_name;
	
	// return created branch
	return new_branch;
}

// remove branch
TreeKinematic::BranchList TreeKinematic::branchRem(const std::string& name) {
	TreeKinematic::BranchList ret_list;
	// trunk cannot be deleted
	if (name.compare(trunk()) == 0) {
		throw(std::runtime_error("Trunk cannot be removed from tree without deleting it."));
	}
	// remove children
	auto it = _branches.find(name);
	if (it == _branches.end()) {
		return ret_list;
	}
	for (auto it_ci: it->second->_children) {
		string child_name = it_ci.second.name;
		TreeKinematic::BranchList child_list = branchRem(it_ci.second.name);
		// append removed children to return list
		for (auto it_cb: child_list) {
			ret_list[child_name] = it_cb.second;
		}
	}
	// remove from parent
	string parent_name = _parent_map[name];
	auto parent_branch = _branches[parent_name];
	parent_branch->_children.erase(name);

	// update parent map
	_parent_map.erase(name);

	// update branches maps
	ret_list[name] = it->second;
	_branches.erase(name);

	return ret_list;
}

// add fruit
Fruit* TreeKinematic::fruitIs(const std::string& name, const std::string& branch_name, double s) {
	// check if fruit exists already.
	// if so, check parent branch
	// if parent is correct, simply return the branch
	// else error out
	auto it_b = _fruit_branch_map.find(name);
	if (it_b != _fruit_branch_map.end()) {
		if (branch_name.compare(it_b->second) == 0) {
			return _fruits[name];
		} else {
			throw(std::runtime_error("Fruit exists with different parent."));
		}
	}

	// check if parent branch exists.
	auto it_p = _branches.find(branch_name);
	if (it_p == _branches.end()) {
		throw(std::runtime_error("Parent branch does not exist."));
	}

	// else add fruit
	Fruit* new_fruit = new Fruit(name);
	// - create child branch info 
	FruitInfo info;
	info.name = name;
	info.s = s;
	// - add to parent branch
	auto br = branch(branch_name);
	br->_fruits[name] = info;
	// - add branch to self
	_fruits[name] = new_fruit;
	// - update parent map
	_fruit_branch_map[name] = branch_name;
	
	// return created branch
	return new_fruit;
}

// remove fruit
Fruit* TreeKinematic::fruitRem(const std::string& name) {
	// remove from parent
	auto it_f = _fruits.find(name);
	if (it_f == _fruits.end()) {
		return NULL;
	}
	auto br = branch(_fruit_branch_map[name]);
	br->_fruits.erase(name);

	// update parent map
	_fruit_branch_map.erase(name);

	// update fruits map
	auto fruit = it_f->second;
	_fruits.erase(name);

	// return removed fruit
	return fruit;
}

// get transform in tree frame
void TreeKinematic::transformInTree(Eigen::Affine3d& ret_trans, const std::string& branch_name, double s) const {
	string parent_name = "";
	ret_trans = Affine3d::Identity();

	string br_name_local = branch_name;
	double s_local = s;
	BranchList::const_iterator br_itr;
	BranchList::const_iterator parent_br_itr;
	Vector3d position;
	Matrix3d rotation;

	// check if valid branch
	br_itr = _branches.find(branch_name);
	if (br_itr == _branches.cend()) {
		throw(runtime_error("Unknown branch."));
	}

	// update transform
	do {
		// get branch
		br_itr = _branches.find(br_name_local);
		// get position in spline
		br_itr->second->spline()->splineOrientation(rotation, s_local);
		br_itr->second->spline()->splineLocation(position, s_local);
		ret_trans.translation() = rotation*ret_trans.translation() + position;
		ret_trans.linear() = rotation * ret_trans.linear();

		// find parent
		if (br_name_local.compare(trunk()) != 0) {
			const auto it = _parent_map.find(branch_name);
			parent_name = it->second;
		} else {
			parent_name = "";
		}

		if (!parent_name.empty()) {
			// get parent
			parent_br_itr = _branches.find(parent_name);
			// get info in parent
			ChildBranchInfo info = parent_br_itr->second->childBranchInfo(br_name_local);
			// update rotation and position
			ret_trans.translation() = info.rotation*ret_trans.translation();
			ret_trans.linear() = info.rotation * ret_trans.linear();
			// update s and br
			br_name_local = parent_name;
			s_local = info.s;
		}
	} while (!parent_name.empty());
}

// get transform in world frame
void TreeKinematic::transformInWorld(Eigen::Affine3d& ret_trans, const std::string& branch_name, double s) const {
	// get orientation in tree
	transformInTree(ret_trans, branch_name, s);
	// apply transform to get position in world
	ret_trans = _transform * ret_trans;
}

// get position in tree frame
void TreeKinematic::positionInTree(Vector3d& ret_vec, const std::string& branch_name, double s) const {
	Affine3d trans;
	transformInTree(trans, branch_name, s);
	ret_vec = trans.translation();
}

// get orientation in tree frame
void TreeKinematic::orientationInTree(Matrix3d& ret_mat, const std::string& branch_name, double s) const {
	Affine3d trans;
	transformInTree(trans, branch_name, s);
	ret_mat = trans.rotation();
}

// get position in world frame
void TreeKinematic::positionInWorld(Vector3d& ret_vec, const std::string& branch_name, double s) const {
	Affine3d trans;
	transformInWorld(trans, branch_name, s);
	ret_vec = trans.translation();
}

// get orientation in world frame
void TreeKinematic::orientationInWorld(Matrix3d& ret_mat, const std::string& branch_name, double s) const {
	Affine3d trans;
	transformInWorld(trans, branch_name, s);
	ret_mat = trans.rotation();
}

// get linear jacobian
void TreeKinematic::jacobianLinear(MatrixXd& ret_mat, const std::string& branch_name, double s) const {

}
