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
	_alpha = 0.0;
	_beta = 0.0;
}

// dtor
QuadraticSplineKinematic::~QuadraticSplineKinematic() {
	// nothing to do
}

// compute position on spline
void QuadraticSplineKinematic::splineLocation(Eigen::Vector3d& ret_vector, double s) const {
	// check for correct ranges
	if(s < 0 || s > _length) { throw(std::runtime_error("splineLocation: Passed value of s is out of bounds.")); }
	_splineLocation(ret_vector, s);
}

// compute orientation on spline
void QuadraticSplineKinematic::splineOrientation(Eigen::Matrix3d& ret_matrix, double s) const {
	// check for correct ranges
	if(s < 0 || s > _length) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }
	_splineOrientation(ret_matrix, s);
}

// get the local position of a point given the deformation coordinates
void QuadraticSplineKinematic::deformedLocation(Eigen::Vector3d& ret_vector, double s, double t, double eta) const {
	// check for correct ranges
	if(s < 0 || s > _length) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }
	if(t < 0 || t > _radius) { throw(std::runtime_error("deformedLocation: Passed value of t is out of bounds.")); }
	if(eta < 0 || eta > 2.0*M_PI) { throw(std::runtime_error("deformedLocation: Passed value of eta is out of bounds.")); }

	// compute the frame position on the spline.
	_splineLocation(ret_vector, s);
	Matrix3d spline_rot;
	_splineOrientation(spline_rot, s);
	Vector3d local_frame_pos(0.0, t*cos(eta), t*sin(eta));
	ret_vector += spline_rot*local_frame_pos;
}

// compute linear jacobian of point on spline
void QuadraticSplineKinematic::splineLinearJacobian(Eigen::MatrixXd& ret_matrix, double s) const {
	if(s < 0 || s > _length) { throw(std::runtime_error("splineOrientation: Passed value of s is out of bounds.")); }

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

// compute jacobian of spline projection length, for spring force
// computations
void QuadraticSplineKinematic::splineTipProjectionLengthJacobian(Eigen::MatrixXd& ret_matrix) const {
	// reshape to 1 x 2
	ret_matrix.setZero(1,2);

	// temp variables
	double g = gam();
	double dg_dalp = dgam_dalp();
	double dg_dbeta = dgam_dbeta();
	double sinc2 = pow(sinc(g/2.0), 2);
	double gam_bracket = sinc2 + g*sinc(g/2.0)*dsinc(g/2.0);

	// lp_dalp
	ret_matrix(0,0) = _length/2.0*gam_bracket*dg_dalp;
	// lp_dbeta
	ret_matrix(0,1) = _length/2.0*gam_bracket*dg_dbeta;
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
