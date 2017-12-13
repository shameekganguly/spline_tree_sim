// QuadraticSplineKinematic.h: Kinematic representation and computations of
// the quadratic spline element

#ifndef QUADRATIC_SPLINE_KINEMATIC_H
#define QUADRATIC_SPLINE_KINEMATIC_H

#include <Eigen/Dense>

// co-ordinates of 3d space within the spline's surface using Cartesian representation
struct SplinePointCartesian {
	double s; // length along the spline
	double py; // y co-ordinate in the planar frame at s
	double pz; // z co-ordinate in the planar frame at s
	// ctor
	SplinePointCartesian (double ss=0.0, double spy=0.0, double spz=0.0)
	: s(ss), py(spy), pz(spz)
	{/* Nothing to do */}
};

// co-ordinates of 3d space within the spline's surface using polar representation
struct SplinePointPolar {
	double s; // length along the spline
	double t; // radial distance in the planar frame at s
	double eta; // angle subtended with X axis in planar frame at s
	// ctor
	SplinePointPolar (double ss=0.0, double st=0.0, double seta=0.0)
	: s(ss), t(st), eta(seta)
	{/* Nothing to do */}
};

class QuadraticSplineKinematic {
// public member functions
public:
	// ctor
	QuadraticSplineKinematic();

	// dtor
	virtual ~QuadraticSplineKinematic();

	// get the local position of a point given the deformation coordinates in polar form
	virtual void deformedLocation(Eigen::Vector3d& ret_vector, const SplinePointPolar& point) const;

	// get the local position of a point given the deformation coordinates in Cartesian form
	virtual void deformedLocation(Eigen::Vector3d& ret_vector, const SplinePointCartesian& point) const;

	// compute position on spline
	virtual void splineLocation(Eigen::Vector3d& ret_vector, double s) const;

	// compute orientation on spline
	virtual void splineOrientation(Eigen::Matrix3d& ret_matrix, double s) const;

	// compute the spline projection length
	virtual double splineTipProjectionLength() const;

	// compute the projection length at a given point along the spline
	virtual double splineProjectionLength(double s) const;

	// compute linear jacobian of point on spline
	virtual void splineLinearJacobian(Eigen::MatrixXd& ret_matrix, double s) const;

	// compute partial derivative of the rotation matrix for the frame at a point along the spline
	virtual void splinedRotdq(Eigen::Matrix3d& ret_dRot_dalp, Eigen::Matrix3d& ret_dRot_dbeta, double s) const;

	// compute linear jacobian of point not necessarily on the spline
	virtual void splineLinearJacobian(Eigen::MatrixXd& ret_matrix, const SplinePointCartesian& spline_point) const;

	// compute jacobian of spline projection length, for spring force
	// computations
	virtual void splineTipProjectionLengthJacobian(Eigen::MatrixXd& ret_matrix) const;

	// compute jacobian of the projection length at a given point along the spline, for spring force
	// computations
	virtual void splineProjectionLengthJacobian(Eigen::MatrixXd& ret_matrix, double s) const;

	// compute curvature
	virtual double splineCurvature() const;

	// compute curvature Jacobian
	virtual void splineCurvatureJacobian(Eigen::MatrixXd& ret_matrix) const;

	// compute closest point on spline to a given point represented in the spline's frame
	// NOTE: this treats the surface of the spline as generated by a constant radius sphere
	// swept along the length of the spline.
	virtual void closestPointToPoint(
		SplinePointCartesian& ret_point,
		double& distance_to_surface,
		Eigen::Vector3d& normal,
		const Eigen::Vector3d& test_point
	) const;

	// compute closest point on spline to a given sphere with the center represented in the spline's frame
	// NOTE: this treats the surface of the spline as generated by a constant radius sphere
	// swept along the length of the spline.
	virtual void closestPointToSphere(
		SplinePointCartesian& ret_point,
		double& distance_to_surface,
		Eigen::Vector3d& normal,
		const Eigen::Vector3d& test_sphere_center,
		const double radius
	) const;

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
