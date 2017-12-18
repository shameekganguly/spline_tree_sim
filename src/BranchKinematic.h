// BranchKinematic.h: Class for data encapsulation and kinematic computations
// for a single branch

#ifndef BRANCH_KINEMATIC_H
#define BRANCH_KINEMATIC_H

#include "QuadraticSplineKinematic.h"
#include "QuadraticSplineDynamic.h"
#include "QuadraticSplineVisual.h"
#include "Fruit.h"

#include <Eigen/Core>
#include <map>
#include <string>
#include <stdexcept>

// child branch info struct
struct ChildBranchInfo {
public:
	// name of child. Empty if struct is not instantiated.
	std::string name;

	// s in parent
	double s;

	// rotation in parent
	Eigen::Matrix3d rotation;

public:
	// ctor
	ChildBranchInfo()
	: name(""), s(0.0), rotation(Eigen::Matrix3d::Identity())
	{ /* nothing to do */ }
};

// fruit info struct
struct FruitInfo {
public:
	// name of fruit. Empty if struct is not instantiated.
	std::string name;

	// s in parent
	double s;

public:
	// ctor
	FruitInfo()
	: name(""), s(0.0)
	{ /* nothing to do */ }
};

// class for data encapsulation and kinematic computations
class BranchKinematic {
public:
	// typedefs
	typedef std::map<std::string, ChildBranchInfo> ChildInfoList;
	typedef std::map<std::string, FruitInfo> FruitInfoList;

	// member functions
public:
	// ctor
	BranchKinematic(const std::string& name)
	: _name(name)
	{
		_spline = new QuadraticSplineKinematic();
		_visual = new QuadraticSplineVisual(_spline);
		_dynamic = new QuadraticSplineDynamic(_spline);
	}

	// dtor
	virtual ~BranchKinematic() {
		delete _spline;
		delete _visual;
		delete _dynamic;
	}

	/* ---- Member access/augment ---- */
	// get name
	const std::string& name() const {
		return _name;
	}

	// get child branch info
	ChildBranchInfo childBranchInfo(const std::string& name) const {
		auto it = _children.find(name);
		if (it != _children.end()) {
			return it->second;
		} else {
			return ChildBranchInfo();
		}
	}

	// update child branch info
	void childBranchInfoIs(const std::string& child_name, const ChildBranchInfo& new_info) {
		// updating name is not allowed
		if (new_info.name != child_name) {
			throw(std::runtime_error("Child name cannot be changed."));
		}

		// if child exists, add new info.
		auto it = _children.find(child_name);
		if (it != _children.end()) {
			it->second = new_info;
		}
	}

	// get non-const iterators
	ChildInfoList::iterator childrenItrBegin() {
		return _children.begin();
	}
	ChildInfoList::iterator childrenItrEnd() {
		return _children.end();
	}

	// get const iterators
	ChildInfoList::const_iterator childrenItrBegin() const {
		return _children.cbegin();
	}
	ChildInfoList::const_iterator childrenItrEnd() const {
		return _children.cend();
	}

	// get fruit info
	FruitInfo fruitInfo(const std::string& name) const {
		auto it = _fruits.find(name);
		if (it != _fruits.end()) {
			return it->second;
		} else {
			return FruitInfo();
		}
	}

	// update fruit info
	void fruitInfoIs(const std::string& fruit_name, const FruitInfo& new_info) {
		// updating name is not allowed
		if (new_info.name != fruit_name) {
			throw(std::runtime_error("Fruit name cannot be changed."));
		}

		// if child exists, add new info.
		auto it = _fruits.find(fruit_name);
		if (it != _fruits.end()) {
			it->second = new_info;
		}
	}

	// get non-const iterators
	FruitInfoList::iterator fruitsItrBegin() {
		return _fruits.begin();
	}
	FruitInfoList::iterator fruitsItrEnd() {
		return _fruits.end();
	}

	// get const iterators
	FruitInfoList::const_iterator fruitsItrBegin() const {
		return _fruits.cbegin();
	}
	FruitInfoList::const_iterator fruitsItrEnd() const {
		return _fruits.cend();
	}

	// get spline kinematic
	QuadraticSplineKinematic* spline() {
		return _spline;
	}

	// get spline kinematic const
	const QuadraticSplineKinematic* spline() const {
		return _spline;
	}

	// get spline visual
	QuadraticSplineVisual* splineVisual() {
		return _visual;
	}

	// get spline visual const
	const QuadraticSplineVisual* splineVisual() const {
		return _visual;
	}

	// get spline dynamic
	QuadraticSplineDynamic* splineDynamic() {
		return _dynamic;
	}

	// get spline dynamic const
	const QuadraticSplineDynamic* splineDynamic() const {
		return _dynamic;
	}

	// data members
public:
	// name of the branch
	std::string _name;

	// list of children of this branch
	ChildInfoList _children;

	// list of fruits of this node
	FruitInfoList _fruits;

	// (kinematic) spline associated with this branch
	QuadraticSplineKinematic* _spline;

	// (visual) spline associated with this branch
	QuadraticSplineVisual* _visual;

	// (dynamic) spline associated with this branch
	QuadraticSplineDynamic* _dynamic;
};

#endif //BRANCH_KINEMATIC_H
