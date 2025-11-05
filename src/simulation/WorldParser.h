#ifndef SIMULATION_WORLD_PARSER_H
#define SIMULATION_WORLD_PARSER_H

#include "TreeParser.h"
#include "World.h"

namespace spline_sim {
struct WorldString {
  std::vector<TreeKinematic *> trees;
};

class WorldParser {
public:
  explicit WorldParser(const std::string &world_file);

  std::unique_ptr<World> LoadWorld();

protected:
  // file name
  std::string file_name_;

  // parsed tree description
  WorldString world_;
};
} // namespace spline_sim

#endif