// A utility file for spline contact implementations
#ifndef SPLINE_CONTACT_H
#define SPLINE_CONTACT_H

#include "QuadraticSplineKinematic.h"
#include "TreeKinematic.h"

#include <Eigen/Core>
#include <string>

// contact between haptic sphere cursor and a branch on the tree
struct ContactInfo {
	// branch name
	std::string branch_name;

	// contact point co-ordinates in spline frame
	SplinePointCartesian point_in_branch;

	// outward normal at contact point in the global frame
	Vector3d normal;

	// penetration depth
	double penetration_depth;

	// ctor
	ContactInfo()
	: branch_name(""), penetration_depth(0.0)
	{/* Nothing to do */}

	// ctor with parameters
	ContactInfo(std::string set_branch_name, SplinePointCartesian set_point_in_branch, double set_penetration_depth)
	: branch_name(set_branch_name), point_in_branch(set_point_in_branch), penetration_depth(set_penetration_depth)
	{/* Nothing to do */}

	// print
	void print() {
		cout << "-- contact info --" << endl;
		cout << "branch name " << branch_name << endl;
		cout << "point in branch: s " << point_in_branch.s 
				<< " py " << point_in_branch.py
				<< " pz " << point_in_branch.pz
				<< endl;
		cout << "normal " << normal << endl;
		cout << "penetration depth " << penetration_depth << endl;
	}
};

// implementation of tree collision detection function
void detectCollisionTreeSphere(
	std::vector<ContactInfo>& ret_list,
	const TreeKinematic* tree,
	const Eigen::Vector3d& sphere_center,
	const double radius
) {
	ret_list.clear();
	TreeKinematic::BranchList::const_iterator br_itr;
	Eigen::Affine3d branch_frame;
	SplinePointCartesian ret_point;
	double distance;
	Vector3d normal;
	for (br_itr=tree->branchesItrBegin(); br_itr!=tree->branchesItrEnd(); ++br_itr) {
		// get position of sphere_center in local branch co-ordinates
		// TODO: this can be made more efficient by getting the branch
		// transforms only locally
		tree->transformInWorld(branch_frame, br_itr->first, 0.0);
		br_itr->second->spline()->closestPointToSphere(
			ret_point,
			distance,
			normal,
			branch_frame.inverse()*sphere_center,
			radius
		);
		// check for collision
		if (distance < 1e-3) {
			ContactInfo info;
			info.branch_name = br_itr->first;
			info.point_in_branch = ret_point;
			info.normal = branch_frame.rotation()*normal;
			info.penetration_depth = (distance < 0)? -distance:0.0;
			ret_list.push_back(info);
		}
	}
}

// struct to encapsulate response of getClosestBranchToPoint
struct CursorDistanceInfo {
	// branch name
	std::string branch_name;

	// contact point co-ordinates in spline frame
	SplinePointCartesian point_in_branch;

	// distance
	double distance;

	// ctor
	CursorDistanceInfo()
	: branch_name(""), distance(0.0)
	{/* Nothing to do */}

	// ctor with parameters
	CursorDistanceInfo(std::string set_branch_name, SplinePointCartesian set_point_in_branch, double set_distance)
	: branch_name(set_branch_name), point_in_branch(set_point_in_branch), distance(set_distance)
	{/* Nothing to do */}

	// print
	void print() {
		cout << "-- cursor distance info --" << endl;
		cout << "branch name " << branch_name << endl;
		cout << "point in branch: s " << point_in_branch.s
				<< " py " << point_in_branch.py
				<< " pz " << point_in_branch.pz
				<< endl;
		cout << "distance " << distance << endl;
	}
};

CursorDistanceInfo getClosestBranchToPoint(const TreeKinematic* tree, const Eigen::Vector3d& point) {
	CursorDistanceInfo ret_info;
	TreeKinematic::BranchList::const_iterator br_itr;
	Eigen::Affine3d branch_frame;
	SplinePointCartesian ret_point;
	double distance;
	Vector3d normal;
	for (br_itr=tree->branchesItrBegin(); br_itr!=tree->branchesItrEnd(); ++br_itr) {
		// get position of point in local branch co-ordinates
		// TODO: this can be made more efficient by getting the branch
		// transforms only locally
		tree->transformInWorld(branch_frame, br_itr->first, 0.0);
		br_itr->second->spline()->closestPointToPoint(
			ret_point,
			distance,
			normal,
			branch_frame.inverse()*point
		);
		// check for minimum distance
		if (ret_info.branch_name.empty() || distance < ret_info.distance) {
			ret_info.branch_name = br_itr->first;
			ret_info.point_in_branch = ret_point;
			ret_info.distance = distance;
		}
	}
	return ret_info;
}

#endif //SPLINE_CONTACT_H
