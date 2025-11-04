// TreeParser: Class to parse tree defined as an xml file

#ifndef TREE_PARSER_H
#define TREE_PARSER_H

#include "BranchKinematic.h"
#include "TreeKinematic.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tinyxml2.h>
#include <vector>

namespace spline_sim {

struct ParentString {
  std::string name;
  double s;
  Eigen::Quaterniond orientation;
  bool f_orientation_assigned;
};

struct SplineString {
  double length;
  double radius;
};

struct DynamicString {
  double ks;
  bool f_ks_assigned;
  double bs;
  bool f_bs_assigned;
  Eigen::Vector3d home_axis;
};

struct BranchString {
  std::string name;
  ParentString parent;
  SplineString spline;
  DynamicString dynamic;
  bool f_dynamic_assigned;
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
  std::string trunk_name;
  bool f_trunk_name_assigned;
  Eigen::Vector3d position;
  bool f_position_assigned;
  Eigen::Quaterniond orientation;
  bool f_orientation_assigned;
  TrunkString trunk;
  std::map<std::string, BranchString> branches;
  std::vector<FruitString> fruits;
  bool fixed = false;
	bool f_fixed_assigned;
};

class TreeParser {
  // member functions
public:
  // ctor
  explicit TreeParser(const std::string &tree_file);

  // load parsed description to a kinematic tree
  // NOTE: the parser DOES NOT deallocate the memory for the created tree
  // upon deletion
  TreeKinematic *loadDescToTree();

  // dtor
  virtual ~TreeParser() = default;

public:
  // file name
  std::string _file_name;

  // parsed tree description
  TreeString _tree;

  // internal member functions
public:
  Eigen::Vector3d parsePosition(const std::string &postion_string);

  Eigen::Quaterniond parseOrientation(const std::string &orientation_string);

  Eigen::Quaterniond parseEulerXZXDeg(const std::string &orientation_string);

  TrunkString parseTrunk(tinyxml2::XMLElement *trunk_element);

  SplineString parseSpline(tinyxml2::XMLElement *spline_element);

  DynamicString parseDynamic(tinyxml2::XMLElement *dynamic_element);

  BranchString parseBranch(tinyxml2::XMLElement *branch_element);

  FruitString parseFruit(tinyxml2::XMLElement *fruit_element);

  ParentString parseParent(tinyxml2::XMLElement *parent_element);
};

} // namespace spline_sim

#endif // TREE_PARSER_H
