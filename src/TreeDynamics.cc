#include "TreeDynamics.h"

#include <Eigen/Core>

using namespace Eigen;

namespace spline_sim {

void TreeDynamics::Step(double dt_sec,
                        const std::vector<ContactInfo> &contact_list,
                        const VectorXd &gamma_external) {
  // MatrixXd J_cs, N;
  // J_cs.setZero(contact_list.size(), tree->dof());
  // N.setZero(contact_list.size(), 3);

  // MatrixXd J_spline_contact_point;
  // VectorXd F_contact, v_contact;
  // VectorXd cc;
  // MatrixXd CM;

  // // Capital gamma denotes net "external" generalized force, which is not
  // // part of the constraint solver step. So not function of q_dot.
  // VectorXd gamma_fruits(tree->dof()), gamma_springs(tree->dof()),
  //     gamma_tree(tree->dof());

  // // - loop over fruits
  // gamma_fruits.setZero();
  // for (const auto fruit_itr = tree->fruitsItrBegin();
  //      fruit_itr != tree->fruitsItrEnd(); ++fruit_itr) {
  //   const std::string branch_name =
  //   tree->fruitParentBranch(fruit_itr->first); double s_cherry =
  //   tree->branch(branch_name)->fruitInfo(fruit_itr->first).s;

  //   Vector3d F_fruit(0, 0, -9.8 * fruit_itr->second->mass());
  //   if (fruit_haptic_force.has_value() &&
  //       fruit_haptic_force->first == fruit_itr->first) {
  //     F_fruit += fruit_haptic_force->second;
  //   }

  //   // TODO: Don't compute on each step
  //   VectorXd Jv_fruit;
  //   tree->jacobianLinear(Jv_fruit, branch_name,
  //                        SplinePointCartesian(s_cherry, 0.0, 0.0));
  //   gamma_fruits += Jv_fruit.transpose() * F_fruit;
  // }
  // // cout << "Gamma cherry" << gamma_fruits.transpose() << endl;

  // // - loop over all branches
  // gamma_springs.setZero();
  // for (branch_itr = tree->branchesItrBegin();
  //      branch_itr != tree->branchesItrEnd(); ++branch_itr) {
  //   branch_ptr = branch_itr->second;
  //   spline_dyn_ptr = branch_ptr->splineDynamic();
  //   // get index in gamma
  //   uint branch_index = tree->branchIndex(branch_ptr->_name);
  //   // compute spring force
  //   spline_dyn_ptr->springForce(branch_spring_force);
  //   gamma_springs.segment<2>(2 * branch_index) = branch_spring_force;

  //   // update tree b_vec
  //   tree_inv_b_vec.segment<2>(2 * branch_index) =
  //       (1.0 / spline_dyn_ptr->_bs) * Vector2d::Ones();
  // }
  // tree_inv_b_mat.diagonal() = tree_inv_b_vec;

  // // - sum fruit and branch spring forces
  // gamma_tree = gamma_fruits + gamma_springs;

  // // compute haptic force on branch
  // if (pulling_branch) {
  //   tree->jacobianLinear(J_branch_pull, closest_cursor_branch_name,
  //                        SplinePointCartesian(closest_cursor_point.s, 0, 0));
  //   gamma_tree += J_branch_pull.transpose() * F_haptic;
  // }

  // // - resolve quasi-static contact
  // if (fHapticDeviceEnabled) {
  //   detectCollisionTreeSphere(contact_list, tree,
  //   cursor->getLocalPos().eigen(),
  //                             2.0 * cursor->getRadius());
  // }

  // if (contact_list.empty()) {
  //   F_proxy_contact.setZero();
  // } else {
  //   // solve simultaneously for the velocities of the splines and the cursor
  //   // - size the contact matrices
  //   J_cs.setZero(contact_list.size(), tree->dof());
  //   N.setZero(contact_list.size(), 3);

  //   // - form the relative velocity Jacobian at the contact points
  //   // NOTE: positive relative velocity along the normal indicates separation
  //   for (uint ind = 0; ind < contact_list.size(); ++ind) {
  //     auto contact_info = contact_list[ind];
  //     // contact_info.print();
  //     // get linear velocity Jacobian from spline
  //     tree->jacobianLinear(J_spline_contact_point, contact_info.branch_name,
  //                          contact_info.point_in_branch);
  //     J_cs.row(ind) << -contact_info.normal.transpose() *
  //                          J_spline_contact_point;
  //     N.row(ind) << contact_info.normal.transpose();
  //   }

  //   // - get the projected dynamics co-efficients
  //   // v_contact = cc + CM*F_contact
  //   cc =
  //       N * (F_proxy - F_haptic) / proxy_b + tree_inv_b_mat * J_cs *
  //       gamma_tree;
  //   CM = N * N.transpose() / proxy_b + tree_inv_b_mat * J_cs *
  //   J_cs.transpose();

  //   // cout << "cc: " << cc.transpose() << endl;
  //   // cout << "CM: " << endl;
  //   // cout << CM << endl;
  //   // - solve the LCP in the contact co-ordinates with the dynamics above
  //   contact_solver.solve(v_contact, F_contact, CM, cc);
  //   // cout << "Solved LCP: " << v_contact.transpose() << " " <<
  //   // F_contact.transpose() << endl;
  //   F_proxy_contact = N.transpose() * F_contact;
  //   gamma_tree += J_cs.transpose() * F_contact;
  // }

  // // - update splines branches
  // for (branch_itr = tree->branchesItrBegin();
  //      branch_itr != tree->branchesItrEnd(); ++branch_itr) {
  //   branch_ptr = branch_itr->second;
  //   spline_ptr = branch_ptr->spline();
  //   spline_dyn_ptr = branch_ptr->splineDynamic();
  //   // get index in gamma
  //   uint branch_index = tree->branchIndex(branch_ptr->_name);
  //   dq[0] = gamma_tree[2 * branch_index + 0] / spline_dyn_ptr->_bs;
  //   dq[1] = gamma_tree[2 * branch_index + 1] / spline_dyn_ptr->_bs;
  //   spline_ptr->_alpha += dq[0] * loop_dt;
  //   spline_ptr->_beta += dq[1] * loop_dt;
  // }
}

void TreeDynamics::Step(double dt_sec,
                        const std::vector<ContactInfo> &contact_list) {
  Step(dt_sec, contact_list, VectorXd::Zero(tree_kinematic_->dof()));
}

void TreeDynamics::Step(double dt_sec,
                        const std::vector<ContactInfo> &contact_list,
                        const FruitHapticForce &fruit_haptic_force) {
  // TODO: compute external Gamma contribution from fruit_haptic_force and call
  // overload with external gamma.
  VectorXd gamma_haptic;
  gamma_haptic.setZero(tree_kinematic_->dof());
  Step(dt_sec, contact_list, gamma_haptic);
}

} // namespace spline_sim