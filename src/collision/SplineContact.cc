#include "SplineContact.h"

#include <iostream>

namespace spline_sim {

void ContactInfo::print() {
  std::cout << "-- contact info --" << "\n";
  std::cout << "branch name " << branch_name << "\n";
  std::cout << "point in branch: s " << point_in_branch.s << " py "
            << point_in_branch.py << " pz " << point_in_branch.pz << "\n";
  std::cout << "normal " << normal << "\n";
  std::cout << "penetration depth " << penetration_depth << std::endl;
}

} // namespace spline_sim