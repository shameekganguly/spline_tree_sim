// QuadraticSplineKinematic.h: Kinematic representation and computations of
// the quadratic spline element

#ifndef QUADRATIC_SPLINE_KINEMATIC_H
#define QUADRATIC_SPLINE_KINEMATIC_H

#include <Eigen/Core>

class QuadraticSplineKinematic {
// public member functions
public:
	// ctor
	QuadraticSplineKinematic();

	// dtor
	virtual ~QuadraticSplineKinematic();

	// get the local position of a point given the deformation coordinates
	virtual void deformedLocation(Eigen::Vector3d& ret_vector, double s, double t, double eta) const;

	// compute position on spline
	virtual void splineLocation(Eigen::Vector3d& ret_vector, double s) const;

	// compute orientation on spline
	virtual void splineOrientation(Eigen::Matrix3d& ret_matrix, double s) const;

	// compute linear jacobian of point on spline
	virtual void splineLinearJacobian(Eigen::MatrixXd& ret_matrix, double s) const;

	// compute jacobian of spline projection length, for spring force
	// computations
	virtual void splineTipProjectionLengthJacobian(Eigen::MatrixXd& ret_matrix) const;

// protected member functions
protected:

// private internal functions
private:
	// compute position on spline. No checks for bounds.
	virtual void _splineLocation(Eigen::Vector3d& ret_vector, double s) const;

	// compute orientation on spline. No checks for bounds.
	virtual void _splineOrientation(Eigen::Matrix3d& ret_matrix, double s) const;

	/* ---- math functions ---- */
	// get gamma
	double gam() const;

	// get derivative of gamma wrt alpha
	double dgam_dalp() const;

	// get derivative of gamma wrt beta
	double dgam_dbeta() const;

	// sinc
	double sinc(double x, double tol=1e-5) const;

	// derivative of sinc
	double dsinc(double x, double tol=1e-5) const;

// data members
// TODO: need to make these protected and add access and write functions, since
// we need to have these parameters be in bounded ranges.
public:
	/* ---- Kinematic info ----*/
	// length of the element
	double _length;

	// radius of the element
	double _radius;

	// alpha angle of the element
	double _alpha;

	// beta angle of the element
	double _beta;
};

#endif //QUADRATIC_SPLINE_KINEMATIC_H
