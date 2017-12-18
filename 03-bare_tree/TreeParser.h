// TreeParser: Class to parse tree defined as an xml file

#ifndef TREE_PARSER_H
#define TREE_PARSER_H

#include "BranchKinematic.h"
#include "TreeKinematic.h"

#include <tinyxml2.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <string>
#include <set>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

struct ParentString {
	std::string name;
	double s;
	Eigen::Quaterniond orientation; bool f_orientation_assigned;
};

struct SplineString {
	double length;
	double radius;
};

struct DynamicString {
	double ks; bool f_ks_assigned;
	double bs; bool f_bs_assigned;
	Eigen::Vector3d home_axis;
};

struct BranchString {
	std::string name;
	ParentString parent;
	SplineString spline;
	DynamicString dynamic; bool f_dynamic_assigned;
};

struct FruitString {
	std::string name;
	double radius;
	ParentString parent;
};

struct TrunkString {
	SplineString spline;
};

struct TreeString {
	std::string name;
	std::string trunk_name; bool f_trunk_name_assigned;
	Eigen::Vector3d position; bool f_position_assigned;
	Eigen::Quaterniond orientation; bool f_orientation_assigned;
	TrunkString trunk;
	std::map<std::string, BranchString> branches;
	std::vector<FruitString> fruits;
};

class TreeParser {
	//member functions
public:
	// ctor
	TreeParser(const std::string& tree_file)
	: _file_name(tree_file)
	{
		std::ifstream model_file (_file_name);

		// reserve memory for the contents of the file
		std::string model_xml_string;
		model_file.seekg(0, std::ios::end);
		model_xml_string.reserve(model_file.tellg());
		model_file.seekg(0, std::ios::beg);
		model_xml_string.assign((std::istreambuf_iterator<char>(model_file)), std::istreambuf_iterator<char>());
		model_file.close();

		tinyxml2::XMLDocument xml_doc;
		xml_doc.Parse(model_xml_string.c_str());
		if (xml_doc.Error()) { throw(std::runtime_error("Couldn't parse xml!")); }

		// parse description to structs
		// - parse tree struct
		tinyxml2::XMLElement *tree_xml = xml_doc.FirstChildElement("tree");
		if(NULL == tree_xml) { throw(std::runtime_error("No tree description found.")); }
		// - - parse tree name
		{
			const char *name = tree_xml->Attribute("name");
			if (!name) { throw(std::runtime_error("Tree must have a name attribute.")); }
			_tree.name = std::string(name);
		}
		// - - parse (optional) trunk name
		{
			const char *name = tree_xml->Attribute("trunk_name");
			if (!name) {
				_tree.f_trunk_name_assigned = false;
			} else {
				_tree.f_trunk_name_assigned = true;
				_tree.trunk_name = name;
			}
		}
		// - - parse optional position in the world
		{
			const char *pos_string = tree_xml->Attribute("position");
			if (!pos_string) {
				_tree.f_position_assigned = false;
			} else {
				_tree.f_position_assigned = true;
				_tree.position = parsePosition(pos_string);
			}
		}
		// - - parse optional orientation in the world
		{
			const char *ori_string = tree_xml->Attribute("orientation");
			if (!ori_string) {
				_tree.f_orientation_assigned = false;
			} else {
				_tree.f_orientation_assigned = true;
				_tree.orientation = parseOrientation(ori_string);
			}
		}

		// - parse trunk struct
		tinyxml2::XMLElement *trunk_xml = tree_xml->FirstChildElement("trunk");
		if(NULL == trunk_xml) { throw(std::runtime_error("No trunk description found.")); }
		_tree.trunk = parseTrunk(trunk_xml);

		// - parse branches
		for (
			tinyxml2::XMLElement* branch_xml = tree_xml->FirstChildElement("branch");
			branch_xml;
			branch_xml = branch_xml->NextSiblingElement("branch")
		) {
			BranchString branch = parseBranch(branch_xml);
			// check for name conflict with existing branches
			auto itr = _tree.branches.find(branch.name);
			if (itr != _tree.branches.end()) {
				throw(std::runtime_error("Branches must have unique names."));
			}
			_tree.branches[branch.name] = branch;
		}

		// - parse fruits
		for (
			tinyxml2::XMLElement* fruit_xml = tree_xml->FirstChildElement("fruit");
			fruit_xml;
			fruit_xml = fruit_xml->NextSiblingElement("fruit")
		) {
			_tree.fruits.push_back(parseFruit(fruit_xml));
		}
	}

	// load parsed description to a kinematic tree
	// NOTE: the parser DOES NOT deallocate the memory for the created tree
	// upon deletion
	TreeKinematic* loadDescToTree() {
		// initialize tree
		TreeKinematic* ret_tree;
		if (_tree.f_trunk_name_assigned) {
			ret_tree = new TreeKinematic(_tree.name, _tree.trunk_name);
		} else {
			ret_tree = new TreeKinematic(_tree.name);
		}
		// add transformation description
		Eigen::Affine3d transform = ret_tree->transformFromWorld();
		if (_tree.f_position_assigned) {
			transform.translation() = _tree.position;
		}
		if (_tree.f_orientation_assigned) {
			transform.linear() = _tree.orientation.toRotationMatrix();
		}
		ret_tree->transformFromWorldIs(transform);
		// add trunk description
		BranchKinematic* trunk_ptr = ret_tree->branch(ret_tree->trunk());
		trunk_ptr->spline()->_length = _tree.trunk.spline.length;
		trunk_ptr->spline()->_radius = _tree.trunk.spline.radius;
		// add branches
		std::set<std::string> branches_added;
		std::map<std::string, BranchString>::iterator bmap_itr;
		bmap_itr = _tree.branches.begin();
		while (branches_added.size() != _tree.branches.size()) {
			if (bmap_itr == _tree.branches.end()) {
				// loop to front
				bmap_itr = _tree.branches.begin();
			}
			// skip if added already
			if (branches_added.count(bmap_itr->first) > 0) {
				++bmap_itr;
			}
			// check if parent is already present
			BranchKinematic* parent = ret_tree->branch(bmap_itr->second.parent.name);
			if (NULL != parent) {
				// insert branch
				auto branch_ptr = ret_tree->branchIs(bmap_itr->second.name, bmap_itr->second.parent.name, 0.0);
				// update branch info in parent
				ChildBranchInfo info;
				info.name = branch_ptr->_name;
				info.s = bmap_itr->second.parent.s;
				if (bmap_itr->second.parent.f_orientation_assigned) {
					info.rotation = bmap_itr->second.parent.orientation.toRotationMatrix();
				}
				ret_tree->branchInfoInParentUpdate(branch_ptr->_name, info);
				// update branch spline info
				branch_ptr->spline()->_length = bmap_itr->second.spline.length;
				branch_ptr->spline()->_radius = bmap_itr->second.spline.radius;
				// update branch spline dynamic info if specified
				if (bmap_itr->second.f_dynamic_assigned) {
					DynamicString dyn_str = bmap_itr->second.dynamic;
					branch_ptr->splineDynamic()->_home_axis = dyn_str.home_axis;
					if (dyn_str.f_ks_assigned) { branch_ptr->splineDynamic()->_ks = dyn_str.ks; }
					if (dyn_str.f_bs_assigned) { branch_ptr->splineDynamic()->_bs = dyn_str.bs; }
					// TODO: also update the home position for the kinematic spline
					// TODO: consider a separate tag for the home position for the kinematic spline
				}
				// move iterator and add to added branches
				branches_added.insert(branch_ptr->_name);
				++bmap_itr;
			}
		}
		// add fruits
		for (auto fdesc_itr: _tree.fruits) {
			// create fruit
			auto fruit = ret_tree->fruitIs(fdesc_itr.name, fdesc_itr.parent.name, 0.0);
			// update fruit info
			fruit->radiusIs(fdesc_itr.radius);
			FruitInfo info;
			info.name = fruit->_name;
			info.s = fdesc_itr.parent.s;
			ret_tree->fruitInfoInParentBranchUpdate(fruit->_name, info);
		}
		// return tree
		return ret_tree;
	}

	// dtor
	virtual ~TreeParser() {
		/* Nothing to do */
	}

public:
	// file name
	std::string _file_name;

	// parsed tree description
	TreeString _tree;

	// internal member functions
public:
	Eigen::Vector3d parsePosition (const std::string& postion_string) {
		Eigen::Vector3d ret_vec;
		std::stringstream ss(postion_string);
		std::string item;
		uint i = 0;
		while (std::getline(ss, item, ' ')) {
			if (i > 2) { throw(std::runtime_error("Malformed position string.")); }
			ret_vec[i++] = std::stod(item);
		}
		return ret_vec;
	}

	Eigen::Quaterniond parseOrientation (const std::string& orientation_string) {
		Eigen::Vector4d temp_vec;
		std::stringstream ss(orientation_string);
		std::string item;
		uint i = 0;
		while (std::getline(ss, item, ' ')) {
			if (i > 3) { throw(std::runtime_error("Malformed orientation string.")); }
			temp_vec[i++] = std::stod(item);
		}
		if(fabs(1.0 - temp_vec.norm()) > 1e-3) { throw(std::runtime_error("Quaternion is not unit norm.")); }
		temp_vec.normalize(); // to remove round off errors
		Eigen::Quaterniond ret_qtn(temp_vec[0], temp_vec[1], temp_vec[2], temp_vec[3]);
		return ret_qtn;
	}

	Eigen::Quaterniond parseEulerXZXDeg (const std::string& orientation_string) {
		Eigen::Vector3d temp_vec;
		std::stringstream ss(orientation_string);
		std::string item;
		uint i = 0;
		while (std::getline(ss, item, ' ')) {
			if (i > 2) { throw(std::runtime_error("Malformed Euler XZX string.")); }
			temp_vec[i++] = std::stod(item);
		}
		Eigen::Quaterniond ret_qtn = Eigen::AngleAxisd(temp_vec[0]*M_PI/180.0, Eigen::Vector3d::UnitX()) *
										Eigen::AngleAxisd(temp_vec[1]*M_PI/180.0, Eigen::Vector3d::UnitZ()) *
										Eigen::AngleAxisd(temp_vec[2]*M_PI/180.0, Eigen::Vector3d::UnitX());
		return ret_qtn;
	}

	TrunkString parseTrunk (tinyxml2::XMLElement* trunk_element) {
		TrunkString ret_trunk;
		tinyxml2::XMLElement *spline_elem = trunk_element->FirstChildElement("spline");
		if(NULL == spline_elem) { throw(std::runtime_error("Spline not specified.")); }
		ret_trunk.spline = parseSpline(spline_elem);
		return ret_trunk;
	}

	SplineString parseSpline (tinyxml2::XMLElement* spline_element) {
		SplineString ret_spline;
		// parse length
		{
			const char *length_str = spline_element->Attribute("length");
			if(!length_str) { throw(std::runtime_error("Missing length.")); }
			ret_spline.length = std::stod(length_str);
		}
		// parse radius
		{
			const char *radius_str = spline_element->Attribute("radius");
			if(!radius_str) { throw(std::runtime_error("Missing radius.")); }
			ret_spline.radius = std::stod(radius_str);
		}
		return ret_spline;
	}

	DynamicString parseDynamic (tinyxml2::XMLElement* dynamic_element) {
		DynamicString ret_dynamic;
		// parse home axis from euler xzx deg. note that the last value is ignored
		// as we only parse the X axis of the resulting frame
		{
			const char *home_str = dynamic_element->Attribute("xzx_deg");
			if(!home_str) { throw(std::runtime_error("Missing home axis definition.")); }
			Eigen::Quaterniond qtn = parseEulerXZXDeg (home_str);
			ret_dynamic.home_axis = qtn.toRotationMatrix().col(0);
		}
		// parse (optional) stiffness parameter ks
		{
			const char *ks_str = dynamic_element->Attribute("ks");
			if(!ks_str) {
				ret_dynamic.f_ks_assigned = false;
			} else {
				ret_dynamic.f_ks_assigned = true;
				ret_dynamic.ks = std::stod(ks_str);
			}
		}
		// parse (optional) damping parameter bs
		{
			const char *bs_str = dynamic_element->Attribute("bs");
			if(!bs_str) {
				ret_dynamic.f_bs_assigned = false;
			} else {
				ret_dynamic.f_bs_assigned = true;
				ret_dynamic.bs = std::stod(bs_str);
			}
		}
		return ret_dynamic;
	}

	BranchString parseBranch (tinyxml2::XMLElement* branch_element) {
		BranchString ret_branch;

		// parse name
		{
			const char *name = branch_element->Attribute("name");
			if (!name) { throw(std::runtime_error("Branch must have a name attribute.")); }
			ret_branch.name = std::string(name);
		}

		// parse parent
		tinyxml2::XMLElement *parent_elem = branch_element->FirstChildElement("parent");
		if(NULL == parent_elem) { throw(std::runtime_error("Branch parent not specified.")); }
		ret_branch.parent = parseParent(parent_elem);

		// parse spline
		tinyxml2::XMLElement *spline_elem = branch_element->FirstChildElement("spline");
		if(NULL == spline_elem) { throw(std::runtime_error("Spline not specified.")); }
		ret_branch.spline = parseSpline(spline_elem);

		// parse optional dynamic
		tinyxml2::XMLElement *dynamic_elem = branch_element->FirstChildElement("dynamic");
		if(NULL == dynamic_elem) {
			ret_branch.f_dynamic_assigned = false;
		} else {
			ret_branch.f_dynamic_assigned = true;
			ret_branch.dynamic = parseDynamic(dynamic_elem);
		}
		ret_branch.spline = parseSpline(spline_elem);

		return ret_branch;
	}

	FruitString parseFruit (tinyxml2::XMLElement* fruit_element) {
		FruitString ret_fruit;

		// parse name
		{
			const char *name = fruit_element->Attribute("name");
			if (!name) { throw(std::runtime_error("Fruit must have a name attribute.")); }
			ret_fruit.name = std::string(name);
		}
		// parse radius
		{
			const char *radius_str = fruit_element->Attribute("radius");
			if (!radius_str) { throw(std::runtime_error("Fruit must have a radius attribute.")); }
			ret_fruit.radius = std::stod(radius_str);
		}
		// parse parent
		tinyxml2::XMLElement *parent_elem = fruit_element->FirstChildElement("parent");
		if(NULL == parent_elem) { throw(std::runtime_error("Fruit parent not specified.")); }
		ret_fruit.parent = parseParent(parent_elem);

		return ret_fruit;
	}

	ParentString parseParent (tinyxml2::XMLElement* parent_element) {
		ParentString ret_parent;
		// parse name
		{
			const char *name = parent_element->Attribute("name");
			if (!name) { throw(std::runtime_error("Parent must have a name attribute.")); }
			ret_parent.name = std::string(name);
		}
		// parse s
		{
			const char *s_string = parent_element->Attribute("s");
			if (!s_string) { throw(std::runtime_error("Parent must have a -s- attribute.")); }
			ret_parent.s = std::stod(s_string);
		}
		// parse optional orientation in quaternion
		bool f_quaternion_specified = false;
		{
			const char *ori_string = parent_element->Attribute("orientation");
			if (!ori_string) {
				ret_parent.f_orientation_assigned = false;
			} else {
				ret_parent.f_orientation_assigned = true;
				f_quaternion_specified = true;
				ret_parent.orientation = parseOrientation(ori_string);
			}
		}
		// parse optional orientation in euler xzx.
		{
			const char *ori_string = parent_element->Attribute("xzx_deg");
			if (!ori_string && !f_quaternion_specified) {
				ret_parent.f_orientation_assigned = false;
			} else if (ori_string && f_quaternion_specified) {
				throw(std::runtime_error("Both quaternion and euler angles specified."));
			}
			else if (ori_string && !f_quaternion_specified) {
				ret_parent.f_orientation_assigned = true;
				ret_parent.orientation = parseEulerXZXDeg(ori_string);
			}
		}

		return ret_parent;
	}
};

#endif //TREE_PARSER_H
