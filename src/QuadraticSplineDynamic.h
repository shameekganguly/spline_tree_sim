// QuadraticSplineDynamic: Dynamics implementation for the spline
#ifndef QUADRATIC_SPLINE_DYNAMIC_H
#define QUADRATIC_SPLINE_DYNAMIC_H

#include "QuadraticSplineKinematic.h"

#include <cmath>

class QuadraticSplineDynamic {
	// member functions
public:
	// ctor
	QuadraticSplineDynamic(const QuadraticSplineKinematic* kinematic)
	: _kinematic(kinematic), _ks(6.0*1e6), _bs(0.50), _home_axis(1.0, 0.0, 0.0) {
		// nothing to do
	}

	// dtor
	virtual ~QuadraticSplineDynamic() {
		// nothing to do
	}

	// get spline stiffness. Note that this might change with alpha and beta
	double stiffness() const {
		return _ks*pow(_kinematic->_radius, 4);
	}

	// get the spring force. the generalized co-ordinates are in the same order
	// as the vector from _kinematic->qVec
	void springForce(Eigen::VectorXd& ret_vec) {
		// // -- curvature based spring force --
		// Eigen::MatrixXd spring_jac;
		// _kinematic->splineCurvatureJacobian(spring_jac);
		// ret_vec = -stiffness()*spring_jac.transpose()*_kinematic->splineCurvature();

		// -- constraint line based spring force --
		Eigen::MatrixXd spring_jac;
		spring_jac.resize(1, 2);
		double sa = sin(_kinematic->_alpha);
		double ca = cos(_kinematic->_alpha);
		double sb = sin(_kinematic->_beta);
		double cb = cos(_kinematic->_beta);
		spring_jac(0,0) = cb*sa*_home_axis[0] + ca*cb*_home_axis[2];
		spring_jac(0,1) = ca*sb*_home_axis[0] - cb*_home_axis[1] - sa*sb*_home_axis[2];
		ret_vec = -stiffness()*spring_jac.transpose()*(1.0 - ca*cb*_home_axis[0] - sb*_home_axis[1] - sa*sb*_home_axis[2]);
	}

	// data members
public:
	// kinematic spline
	const QuadraticSplineKinematic* _kinematic;

	// stiffness constants
	double _ks;

	// damping constants (for first order dynamics)
	double _bs;

	// home axis specification
	Eigen::Vector3d _home_axis;

	// TODO: inertial dynamics

	// TODO: gravitational dynamics
};

#endif //QUADRATIC_SPLINE_DYNAMIC_H
