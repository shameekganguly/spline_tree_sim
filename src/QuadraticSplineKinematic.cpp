#include "QuadraticSplineKinematic.h"
#include <math.h>
#include <stdexcept>

#include <iostream>
using namespace std;

using namespace Eigen;

// ctor
QuadraticSplineKinematic::QuadraticSplineKinematic()
{
	// default values
	_radius = 0.04;
	_length = 1.0;
	_alpha = 0.00001;
	_beta = 0.00001;
}

// dtor
QuadraticSplineKinematic::~QuadraticSplineKinematic() {
	// nothing to do
}

// compute position on spline
void QuadraticSplineKinematic::splineLocation(Eigen::Vector3d& ret_vector, double s) const {
	// check for correct ranges
	if(s < (0 - 1e-10) || s > (_length+1e-10)) { throw(std::runtime_error("splineLocation: Passed value of s is out of bounds.")); }
	_splineLocation(ret_vector, s);
}

// compute orientation on spline
void QuadraticSplineKinematic::splineOrientation(Eigen::Matrix3d& ret_matrix, double s) const {
	// check for correct ranges
	if(s < (0 - 1e-10) || s > (_length+1e-10)) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }
	_splineOrientation(ret_matrix, s);
}

// get the local position of a point given the deformation coordinates in polar form
void QuadraticSplineKinematic::deformedLocation(Eigen::Vector3d& ret_vector, const SplinePointPolar& point) const {
	// check for correct ranges
	if(point.s < (0 - 1e-10) || point.s > (_length+1e-10)) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }
	if(point.t < (0 - 1e-10) || point.t > (_radius+1e-10)) { throw(std::runtime_error("deformedLocation: Passed value of t is out of bounds.")); }
	if(point.eta < (0 - 1e-10) || point.eta > (2.0*M_PI+1e-10)) { throw(std::runtime_error("deformedLocation: Passed value of eta is out of bounds.")); }

	// compute the frame position on the spline.
	_splineLocation(ret_vector, point.s);
	Matrix3d spline_rot;
	_splineOrientation(spline_rot, point.s);
	Vector3d local_frame_pos(0.0, point.t*cos(point.eta), point.t*sin(point.eta));
	ret_vector += spline_rot*local_frame_pos;
}

// get the local position of a point given the deformation coordinates in Cartesian form
void QuadraticSplineKinematic::deformedLocation(Eigen::Vector3d& ret_vector, const SplinePointCartesian& point) const {
	double t = sqrt(point.py*point.py + point.pz*point.pz);
	// check for correct ranges
	if(point.s < (0 - 1e-10) || point.s > (_length+1e-10)) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }
	if(t < (0 - 1e-10) || t > (_radius+1e-10)) { throw(std::runtime_error("deformedLocation: Passed value of px and py are out of bounds.")); }

	// compute the frame position on the spline.
	_splineLocation(ret_vector, point.s);
	Matrix3d spline_rot;
	_splineOrientation(spline_rot, point.s);
	Vector3d local_frame_pos(0.0, point.py, point.pz);
	ret_vector += spline_rot*local_frame_pos;
}

// compute the spline projection length
double QuadraticSplineKinematic::splineTipProjectionLength() const {
	return splineProjectionLength(_length);
}

// compute the projection length at a given point along the spline
double QuadraticSplineKinematic::splineProjectionLength(double s) const {
	return s*s*gam()/(2.0*_length)*pow(sinc(s*gam()/2.0/_length) ,2);
}

// compute linear jacobian of point on spline
void QuadraticSplineKinematic::splineLinearJacobian(Eigen::MatrixXd& ret_matrix, double s) const {
	if(s < (0 - 1e-10) || s > (_length+1e-10)) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }

	// reshape to 3 x 2
	ret_matrix.setZero(3,2);

	// temp variables
	double g = gam();
	double tb = tan(_beta);
	double sa = sin(_alpha);
	double ca = cos(_alpha);
	double cb = cos(_beta);
	double den = sqrt(tb*tb + sa*sa);
	double dg_dalp = dgam_dalp();
	double dg_dbeta = dgam_dbeta();
	double sinc2 = pow(sinc(s*g/(2.0*_length)), 2);
	double gam_bracket = sinc2 + s*g/_length*sinc(s*g/(2.0*_length))*dsinc(s*g/(2.0*_length));

	// px_dalp
	ret_matrix(0,0) = s*s/_length*dsinc(s*g/_length)*dg_dalp;

	// py_dalp
	ret_matrix(1,0) = s*s/(2.0*_length)*gam_bracket*tb/den*dg_dalp
						- s*s*g/(2.0*_length)*sinc2*tb*sa*ca/pow(den, 3);

	// pz_dalp
	ret_matrix(2,0) = -s*s/(2.0*_length)*gam_bracket*sa/den*dg_dalp
						- s*s*g/(2.0*_length)*sinc2*tb*tb*ca/pow(den, 3);

	// px_dbeta
	ret_matrix(0,1) = s*s/_length*dsinc(s*g/_length)*dg_dbeta;

	// py_dbeta
	ret_matrix(1,1) = s*s/(2.0*_length)*gam_bracket*tb/den*dg_dbeta
						+ s*s*g/(2.0*_length)*sinc2*sa*sa/pow(den, 3)/pow(cb, 2);

	// pz_dbeta
	ret_matrix(2,1) = -s*s/(2.0*_length)*gam_bracket*sa/den*dg_dbeta
						+ s*s*g/(2.0*_length)*sinc2*sa*tb/pow(den, 3)/pow(cb, 2);
}

// compute linear jacobian of point not necessarily on the spline
void QuadraticSplineKinematic::splineLinearJacobian(Eigen::MatrixXd& ret_matrix, const SplinePointCartesian& spline_point) const {
	double s = spline_point.s;
	double py = spline_point.py;
	double pz = spline_point.pz;

	double t = sqrt(py*py + pz*pz);
	// check for correct ranges
	if(s < (0 - 1e-10) || s > (_length+1e-10)) { throw(std::runtime_error("Passed value of s is out of bounds.")); }
	if(t < (0 - 1e-10) || t > (_radius+1e-10)) { throw(std::runtime_error("Passed value of px and py are out of bounds.")); }

	Matrix3d dRot_dalp, dRot_dbeta;
	MatrixXd J_vs;
	Vector3d local_frame_pos(0.0, py, pz);
	splineLinearJacobian(J_vs, s);

	// if point is very close to the spline, simply return this
	if (t < 1e-3) { //TODO: normalize by radius and length
		ret_matrix = J_vs;
		return;
	}
	splinedRotdq(dRot_dalp, dRot_dbeta, s);
	ret_matrix.setZero(3, 2);
	ret_matrix.col(0) = dRot_dalp*local_frame_pos + J_vs.col(0);
	ret_matrix.col(1) = dRot_dbeta*local_frame_pos + J_vs.col(1);
}

// compute partial derivative of the rotation matrix for the frame at a point along the spline
void QuadraticSplineKinematic::splinedRotdq(Eigen::Matrix3d& ret_dRot_dalp, Eigen::Matrix3d& ret_dRot_dbeta, double s) const {
	if(s < (0 - 1e-10) || s > (_length+1e-10)) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }

	// temp variables
	double g = gam();
	double tb = tan(_beta);
	double sa = sin(_alpha);
	double sb = sin(_beta);
	double ca = cos(_alpha);
	double cb = cos(_beta);
	double den = sqrt(tb*tb + sa*sa);
	double dg_dalp = dgam_dalp();
	double dg_dbeta = dgam_dbeta();
	double sinc2 = pow(sinc(s*g/(2.0*_length)), 2);
	double gam_bracket = sinc2 + s*g/_length*sinc(s*g/(2.0*_length))*dsinc(s*g/(2.0*_length));

	/* ---- dRot_dalp ---- */
	// Rxx_dalp
	ret_dRot_dalp(0,0) = -s/_length*sin(s*g/_length)*dg_dalp;
	// Rxy_dalp
	ret_dRot_dalp(1,0) = s/_length*cos(s*g/_length)*tb/den*dg_dalp;
	ret_dRot_dalp(1,0) -= sin(s*g/_length)*(ca*sa*tb)/pow(den, 3);
	// Rxz_dalp
	ret_dRot_dalp(2,0) = -s/_length*cos(s*g/_length)*sa/den*dg_dalp;
	ret_dRot_dalp(2,0) -= sin(s*g/_length)*ca*(tb*tb)/pow(den, 3);
	// Ryx_dalp
	ret_dRot_dalp(0,1) = -ret_dRot_dalp(1,0);
	// Ryy_dalp
	ret_dRot_dalp(1,1) = 2.0*ca*sa*(1.0 - cos(s*g/_length))/pow(den, 2);
	ret_dRot_dalp(1,1) -= s/_length*sin(s*g/_length)*dg_dalp;
	ret_dRot_dalp(1,1) *= tb*tb/pow(den, 2);
	// Ryz_dalp
	ret_dRot_dalp(2,1) = s/_length*sin(s*g/_length)*sa*tb/(tb*tb + sa*sa)*dg_dalp;
	ret_dRot_dalp(2,1) += (1.0 - cos(s*g/_length))*ca*tb*(tb*tb - sa*sa)/pow(den, 4);
	// Rzx_dalp
	ret_dRot_dalp(0,2) = -ret_dRot_dalp(2,0);
	// Rzy_dalp
	ret_dRot_dalp(1,2) = ret_dRot_dalp(2,1);
	// Rzz_dalp
	ret_dRot_dalp(2,2) = 2.0*ca*(tb*tb)/sa*(cos(s*g/_length) - 1)/pow(den, 2);
	ret_dRot_dalp(2,2) -= s/_length*sin(s*g/_length)*dg_dalp;
	ret_dRot_dalp(2,2) *= sa*sa/pow(den, 2);

	/* ---- dRot_dbeta ---- */
	// Rxx_dbeta
	ret_dRot_dbeta(0,0) = -s/_length*sin(s*g/_length)*dg_dbeta;
	// Rxy_dbeta
	ret_dRot_dbeta(1,0) = -s/_length*cos(s*g/_length)*tb/den*dg_dbeta;
	ret_dRot_dbeta(1,0) -= sin(s*g/_length)*(sa*sa)/(cb*cb)/pow(den, 3);
	// Rxz_dbeta
	ret_dRot_dbeta(2,0) = s/_length*cos(s*g/_length)*sa/den*dg_dbeta;
	ret_dRot_dbeta(2,0) -= sin(s*g/_length)*sa*tb/(cb*cb)/pow(den, 3);
	// Ryx_dbeta
	ret_dRot_dbeta(0,1) = -ret_dRot_dbeta(1,0);
	// Ryy_dbeta
	ret_dRot_dbeta(1,1) = 2.0*(sa*sa)/(cb*sb)*(cos(s*g/_length)-1.0)/pow(den, 2);
	ret_dRot_dbeta(1,1) -= s/_length*sin(s*g/_length)*dg_dbeta;
	ret_dRot_dbeta(1,1) *= tb*tb/pow(den, 2);
	// Ryz_dbeta
	ret_dRot_dbeta(2,1) = s/_length*sin(s*g/_length)*sa*tb/pow(den, 2)*dg_dbeta;
	ret_dRot_dbeta(2,1) += (1.0 - cos(s*g/_length))*sa/(cb*cb)*(sa*sa - tb*tb)/pow(den, 4);
	// Rzx_dbeta
	ret_dRot_dbeta(0,2) = -ret_dRot_dbeta(2,0);
	// Rzy_dbeta
	ret_dRot_dbeta(1,2) = ret_dRot_dbeta(2,1);
	// Rzz_dbeta
	ret_dRot_dbeta(2,2) = 2.0*tb/(cb*cb)*(1.0-cos(s*g/_length))/pow(den, 2);
	ret_dRot_dbeta(2,2) -= s/_length*sin(s*g/_length)*dg_dbeta;
	ret_dRot_dbeta(2,2) *= sa*sa/pow(den, 2);
}

// compute jacobian of spline projection length, for spring force
// computations
void QuadraticSplineKinematic::splineTipProjectionLengthJacobian(Eigen::MatrixXd& ret_matrix) const {
	splineProjectionLengthJacobian(ret_matrix, _length);
}

// compute jacobian of the projection length at a given point along the spline, for spring force
// computations
void QuadraticSplineKinematic::splineProjectionLengthJacobian(Eigen::MatrixXd& ret_matrix, double s) const {
	// reshape to 1 x 2
	ret_matrix.setZero(1,2);

	// temp variables
	double g = gam();
	double dg_dalp = dgam_dalp();
	double dg_dbeta = dgam_dbeta();
	double sinc2 = pow(sinc(s*g/(2.0*_length)), 2);
	double gam_bracket = sinc2 + s*g/_length*sinc(s*g/(2.0*_length))*dsinc(s*g/(2.0*_length));

	// lp_dalp
	ret_matrix(0,0) = s*s*g/(2.0*_length)*gam_bracket*dg_dalp;
	// lp_dbeta
	ret_matrix(0,1) = s*s*g/(2.0*_length)*gam_bracket*dg_dbeta;
}

// internal: compute position on spline
void QuadraticSplineKinematic::_splineLocation(Eigen::Vector3d& ret_vector, double s) const {
	ret_vector[0] = s*sinc(s*gam()/_length);
	if (fabs(_alpha) < 1e-3 && fabs(_beta) < 1e-3) {
		ret_vector[1] = 0.0;
		ret_vector[2] = 0.0;
	}
	else {
		double sp = s*s*gam()/(2.0*_length) * pow(sinc(s*gam()/(2.0*_length)), 2);
		double den = sqrt(pow(tan(_beta),2) + pow(sin(_alpha),2));
		ret_vector[1] = sp * tan(_beta)/den;
		ret_vector[2] = sp * -sin(_alpha)/den;
	}
}

// internal: compute orientation on spline
void QuadraticSplineKinematic::_splineOrientation(Eigen::Matrix3d& ret_matrix, double s) const {
	if (fabs(_alpha) < 1e-3 && fabs(_beta) < 1e-3) {
		ret_matrix.setIdentity();
	} else {
		double phi = s*gam()/_length;
		double den = sqrt(pow(tan(_beta),2) + pow(sin(_alpha),2));
		double cos_t = tan(_beta)/den;
		double sin_t = -sin(_alpha)/den;
		ret_matrix << cos(phi), -sin(phi)*cos_t, -sin(phi)*sin_t,
					sin(phi)*cos_t, pow(sin_t,2) + pow(cos_t,2)*cos(phi), (cos(phi)-1.0)*cos_t*sin_t,
					sin(phi)*sin_t, (cos(phi)-1.0)*cos_t*sin_t, pow(cos_t,2) + pow(sin_t,2)*cos(phi);
	}
}

// internal: get gamma for current alpha, beta
double QuadraticSplineKinematic::gam() const {
	return M_PI - 2.0*asin(cos(_alpha)*cos(_beta));
}

// get derivative of gamma wrt alpha
double QuadraticSplineKinematic::dgam_dalp() const {
	double psi = asin(cos(_alpha)*cos(_beta));
	double dpsi_dalp = -cos(_beta)*sin(_alpha)/cos(psi);
	return -2.0*dpsi_dalp;
}

// get derivative of gamma wrt beta
double QuadraticSplineKinematic::dgam_dbeta() const {
	double psi = asin(cos(_alpha)*cos(_beta));
	double dpsi_dbeta = -sin(_beta)*cos(_alpha)/cos(psi);
	return -2.0*dpsi_dbeta;
}

// internal: sinc function
double QuadraticSplineKinematic::sinc(double x, double tol) const {
	if (fabs(x) < tol) { return 1.0; }
	return sin(x)/x;
}

// internal: derivative sinc function
double QuadraticSplineKinematic::dsinc(double x, double tol) const {
	if (fabs(x) < tol) { return 0.0; }
	return (cos(x) - sinc(x))/x;
}
