#include "WorldParser.h"

namespace spline_sim {
WorldParser::WorldParser(const std::string &world_file)
    : file_name_(world_file) {
  std::ifstream model_file(file_name_);

  // reserve memory for the contents of the file
  std::string model_xml_string;
  model_file.seekg(0, std::ios::end);
  model_xml_string.reserve(model_file.tellg());
  model_file.seekg(0, std::ios::beg);
  model_xml_string.assign((std::istreambuf_iterator<char>(model_file)),
                          std::istreambuf_iterator<char>());
  model_file.close();

  tinyxml2::XMLDocument xml_doc;
  xml_doc.Parse(model_xml_string.c_str());
  if (xml_doc.Error()) {
    throw(std::runtime_error("Couldn't parse xml!"));
  }

  tinyxml2::XMLElement *world_xml = xml_doc.FirstChildElement("world");
  if (world_xml == nullptr) {
    throw(std::runtime_error("No world element found."));
  }

  // parse description to structs
  // - parse tree structs
  for (tinyxml2::XMLElement *tree_xml = world_xml->FirstChildElement("tree");
       tree_xml; tree_xml = tree_xml->NextSiblingElement("tree")) {
    // - parse tree from given path
    tinyxml2::XMLElement *path_element = tree_xml->FirstChildElement("path");
    if (!path_element) {
      throw(std::runtime_error("World::Tree must have a path element."));
    }
    TreeParser parser(std::string(path_element->GetText()));
    TreeKinematic *tree_kinematic = parser.loadDescToTree();

    // - parse tree name
    {
      const char *name = tree_xml->Attribute("name");
      if (!name) {
        throw(std::runtime_error("World::Tree must have a name attribute."));
      }
      tree_kinematic->nameIs(std::string(name));
    }

    // - parse tree transform
    // - - parse optional position in the world
    auto world_T_tree = Eigen::Affine3d::Identity();
    {
      const char *pos_string = tree_xml->Attribute("position");
      if (pos_string) {
        world_T_tree.translation() = TreeParser::parsePosition(pos_string);
      }
    }
    // - - parse optional orientation in the world
    {
      const char *ori_string = tree_xml->Attribute("orientation");
      if (ori_string) {
        world_T_tree.linear() =
            TreeParser::parseOrientation(ori_string).toRotationMatrix();
      }
    }
    tree_kinematic->transformFromWorldIs(world_T_tree);

    // - parse optional fixed attribute
    {
      const char *bool_string = tree_xml->Attribute("fixed");
      if (bool_string) {
        tree_kinematic->FixedIs(TreeParser::parseBool(bool_string));
      }
    }
    world_.trees.push_back(tree_kinematic);
  }
  if (world_.trees.empty()) {
    throw(std::runtime_error("World does not have any trees."));
  }
}

std::unique_ptr<World> WorldParser::LoadWorld() {
  auto world = std::make_unique<World>();
  for (auto tree_ptr : world_.trees) {
    world->AddTree(std::unique_ptr<TreeKinematic>(tree_ptr));
  }
  return world;
}
} // namespace spline_sim