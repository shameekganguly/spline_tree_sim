// QuadraticSplineVisual.h: Chai visual for a single quadratic spline element

#ifndef QUADRATIC_SPLINE_VISUAL_H
#define QUADRATIC_SPLINE_VISUAL_H

#include <Eigen/Core>
#include <chai3d.h>

// NOTE:
// Pose of frame with respect to parent: The origin is located at the base of
// the spline, in the center of the plane. The X axis points towards the length
// of the spline. Y and Z axes are default, orthogonal to X.

class QuadraticSplineVisual: public chai3d::cGenericObject {
// public member functions
public:
	// ctor
	QuadraticSplineVisual();

	// dtor
	virtual ~QuadraticSplineVisual();

// protected member functions
protected:
	// render: from parent class
	void render(chai3d::cRenderOptions& a_options);

// data members
// TODO: need to make these protected and add access and write functions, since
// we need to have these parameters be in bounded ranges.
public:
	/* ---- Kinematic info ----*/
	// radius of the element
	double _radius;

	// length of the element
	double _length;

	// alpha angle of the element
	double _alpha;

	// beta angle of the element
	double _beta;

	/* ---- Graphic info ----*/
protected:
	// spacing between vertex sets on plane cross sections along the length of the spline
	double _ds_longitudinal;

	// subtended angle between vertices in a single plane cross section
	double _ds_angle;

public:
	// number of vertices in a single plane cross section
	uint _nv_plane;

	// number of plane cross sections along the length of the spline
	uint _nv_longitudinal;
};

#endif //QUADRATIC_SPLINE_VISUAL_H
