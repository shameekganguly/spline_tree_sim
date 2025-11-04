#include "TreeDynamics.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

namespace spline_sim {

namespace {
constexpr double kFloatingBaseLinearDamping = 0.1;
constexpr double kFloatingBaseAngularDamping = 0.01;
} // namespace

void TreeDynamics::Step(double dt_sec,
                        const std::vector<ContactInfo> &contact_list,
                        const VectorXd &gamma_external) {
  // TODO: make these class-internal to save some runtime memory allocation.
  // Capital gamma denotes net "external" generalized force, which is not
  // part of the constraint solver step. So not function of q_dot.
  VectorXd gamma_fruits(tree_kinematic_->dof()),
      gamma_springs(tree_kinematic_->dof()), gamma_tree(tree_kinematic_->dof());
  VectorXd tree_inv_b_vec(tree_kinematic_->dof());
  const uint base_fixed_col_offset = tree_kinematic_->Fixed() ? 0 : 6;

  // - loop over fruits
  gamma_fruits.setZero();
  for (auto fruit_itr = tree_kinematic_->fruitsItrBegin();
       fruit_itr != tree_kinematic_->fruitsItrEnd(); ++fruit_itr) {
    const std::string branch_name =
        tree_kinematic_->fruitParentBranch(fruit_itr->first);
    double s_cherry =
        tree_kinematic_->branch(branch_name)->fruitInfo(fruit_itr->first).s;

    Vector3d F_fruit(0, 0, -9.8 * fruit_itr->second->mass());

    // TODO: Don't compute on each step
    MatrixXd Jv_fruit;
    tree_kinematic_->jacobianLinear(Jv_fruit, branch_name,
                                    SplinePointCartesian(s_cherry, 0.0, 0.0));
    gamma_fruits += Jv_fruit.transpose() * F_fruit;
  }
  // cout << "Gamma cherry" << gamma_fruits.transpose() << endl;

  // - loop over all branches
  gamma_springs.setZero();
  tree_inv_b_vec.setZero();
  for (auto branch_itr = tree_kinematic_->branchesItrBegin();
       branch_itr != tree_kinematic_->branchesItrEnd(); ++branch_itr) {
    auto *branch_ptr = branch_itr->second;
    auto *spline_dyn_ptr = branch_ptr->splineDynamic();
    // get index in gamma
    uint branch_index = tree_kinematic_->branchIndex(branch_ptr->_name);
    // compute spring force
    Eigen::Vector2d branch_spring_force;
    spline_dyn_ptr->springForce(branch_spring_force);
    gamma_springs.segment<2>(base_fixed_col_offset + 2 * branch_index) =
        branch_spring_force;

    // update tree b_vec
    tree_inv_b_vec.segment<2>(base_fixed_col_offset + 2 * branch_index)
        .fill(1.0 / spline_dyn_ptr->_bs);
  }
  DiagonalMatrix<double, -1> tree_inv_b_mat(tree_kinematic_->dof());
  tree_inv_b_mat.diagonal() = tree_inv_b_vec;

  // - sum fruit and branch spring forces
  gamma_tree = gamma_fruits + gamma_springs + gamma_external;

  // - Add base gravity force if it is not fixed
  if (!tree_kinematic_->Fixed()) {
    auto trunk_spline_dyn_ptr =
        tree_kinematic_->branch(tree_kinematic_->trunk())->splineDynamic();
    gamma_tree[2] = -9.8 * trunk_spline_dyn_ptr->_mass;

    tree_inv_b_mat.diagonal().segment<3>(0).fill(1.0 /
                                                 kFloatingBaseLinearDamping);
    tree_inv_b_mat.diagonal().segment<3>(3).fill(1.0 /
                                                 kFloatingBaseLinearDamping);
  }
  // TODO: Add branch gravity force

  // -- TODO: Move to tree haptic controller
  // // compute haptic force on branch
  // if (pulling_branch) {
  //   tree_kinematic_->jacobianLinear(J_branch_pull,
  //   closest_cursor_branch_name,
  //                        SplinePointCartesian(closest_cursor_point.s, 0, 0));
  //   gamma_tree += J_branch_pull.transpose() * F_haptic;
  // }

  // // - detect collision between tree and haptic cursor
  // if (fHapticDeviceEnabled) {
  //   detectCollisionTreeSphere(contact_list, tree,
  //   cursor->getLocalPos().eigen(),
  //                             2.0 * cursor->getRadius());
  // }
  // --

  if (!contact_list.empty()) {
    // Contact space Jacobian
    MatrixXd J_cs;

    // Normal unit vectors, row-wise
    MatrixXd N;

    // solve simultaneously for the velocities of the splines and the cursor
    // - size the contact matrices
    J_cs.setZero(contact_list.size(), tree_kinematic_->dof());
    N.setZero(contact_list.size(), 3);

    // - form the relative velocity Jacobian at the contact points
    // NOTE: positive relative velocity along the normal indicates separation
    for (uint ind = 0; ind < contact_list.size(); ++ind) {
      auto contact_info = contact_list[ind];
      // contact_info.print();
      // get linear velocity Jacobian from spline
      MatrixXd J_spline_contact_point;
      tree_kinematic_->jacobianLinear(J_spline_contact_point,
                                      contact_info.branch_name,
                                      contact_info.point_in_branch);
      J_cs.row(ind) << -contact_info.normal.transpose() *
                           J_spline_contact_point;
      N.row(ind) << contact_info.normal.transpose();
    }

    // Contact dynamics (without cursor):
    //     Bmat * q_dot = Gamma + J_cs^T * F_contact
    // =>  q_dot = Bmat_inv * Gamma + Bmat_inv * J_cs^T * F_contact
    // =>  v_contact =
    //         J_cs * Bmat_inv * Gamma + J_cs * Bmat_inv * J_cs^T * F_contact
    // define cc = J_cs * Bmat_inv * Gamma  is the un-forced contact velocity
    // define CM = J_cs * Bmat_inv * J_cs^T

    // - get the projected dynamics co-efficients
    // v_contact = cc + CM*F_contact

    VectorXd cc = J_cs * tree_inv_b_mat * gamma_tree;
    MatrixXd CM = J_cs * tree_inv_b_mat * J_cs.transpose();

    // TODO: - add haptic proxy dynamics with callback:
    // VectorXd cc =
    //     N * (F_proxy - F_haptic) / proxy_b + J_cs * tree_inv_b_mat *
    //     gamma_tree;
    // MatrixXd CM =
    //     N * N.transpose() / proxy_b + J_cs * tree_inv_b_mat *
    //     J_cs.transpose();

    // cout << "cc: " << cc.transpose() << endl;
    // cout << "CM: " << endl;
    // cout << CM << endl;
    // - solve the LCP in the contact co-ordinates with the dynamics above
    VectorXd F_contact, v_contact;
    contact_solver.solve(v_contact, F_contact, CM, cc);
    // cout << "Solved LCP: " << v_contact.transpose() << " " <<
    // F_contact.transpose() << endl;
    gamma_tree += J_cs.transpose() * F_contact;

    // TODO: set contact force on haptic proxy with callback.
    // F_proxy_contact = N.transpose() * F_contact;
  }

  // - update splines branches
  for (auto branch_itr = tree_kinematic_->branchesItrBegin();
       branch_itr != tree_kinematic_->branchesItrEnd(); ++branch_itr) {
    auto *branch_ptr = branch_itr->second;
    auto *spline_ptr = branch_ptr->spline();
    auto *spline_dyn_ptr = branch_ptr->splineDynamic();

    // Get index in gamma
    uint branch_index = tree_kinematic_->branchIndex(branch_ptr->_name);

    // Compute gc velocity
    double dalpha_dt, dbeta_dt;
    dalpha_dt = gamma_tree[base_fixed_col_offset + 2 * branch_index] /
                spline_dyn_ptr->_bs;
    dbeta_dt = gamma_tree[base_fixed_col_offset + 2 * branch_index + 1] /
               spline_dyn_ptr->_bs;

    // Integrate
    spline_ptr->_alpha += dalpha_dt * dt_sec;
    spline_ptr->_beta += dbeta_dt * dt_sec;
    // std::cout << branch_itr->first << " a: " << spline_ptr->_alpha
    //           << " b: " << spline_ptr->_beta << " gamma: " <<
    //           spline_ptr->gam()
    //           << "\n";
  }

  if (!tree_kinematic_->Fixed()) {
    // Update world transform
    Vector3d base_v = gamma_tree.segment<3>(0) / kFloatingBaseLinearDamping;
    Vector3d base_w = gamma_tree.segment<3>(3) / kFloatingBaseAngularDamping;

    Affine3d new_transform = tree_kinematic_->transformFromWorld();
    new_transform.pretranslate(base_v * dt_sec);
    if (base_w.norm() > 1e-6) {
      new_transform.prerotate(
          AngleAxisd(base_w.norm() * dt_sec, base_w / base_w.norm()));
    }
    tree_kinematic_->transformFromWorldIs(new_transform);
  }
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

  for (auto fruit_itr = tree_kinematic_->fruitsItrBegin();
       fruit_itr != tree_kinematic_->fruitsItrEnd(); ++fruit_itr) {
    if (fruit_haptic_force.first == fruit_itr->first) {
      const std::string branch_name =
          tree_kinematic_->fruitParentBranch(fruit_itr->first);
      double s_fruit =
          tree_kinematic_->branch(branch_name)->fruitInfo(fruit_itr->first).s;

      // TODO: Don't compute on each step
      MatrixXd Jv_fruit;
      tree_kinematic_->jacobianLinear(Jv_fruit, branch_name,
                                      SplinePointCartesian(s_fruit, 0.0, 0.0));
      gamma_haptic = Jv_fruit.transpose() * fruit_haptic_force.second;
      break;
    }
  }

  Step(dt_sec, contact_list, gamma_haptic);
}

} // namespace spline_sim