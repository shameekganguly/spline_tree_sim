// Fruit.h: Class encapsulating a single fruit's data

#ifndef FRUIT_H
#define FRUIT_H

#include <world/CShapeSphere.h>

#include <string>
#include <cmath>
#include <stdexcept>

class Fruit {
	// member functions
public:
	// ctor
	Fruit(const std::string& name, double radius=0.04, double density=0.65)
	: _name(name), _density(density)
	{
		// create a sphere to represent this fruit
		_graphic = new chai3d::cShapeSphere(radius);
	}

	// dtor
	virtual ~Fruit() {
		delete _graphic;
	}

	// get graphic
	chai3d::cShapeSphere* graphic() {
		return _graphic;
	}

	// get graphic const
	const chai3d::cShapeSphere* graphic() const {
		return _graphic;
	}

	// get radius
	double radius() const {
		return _graphic->getRadius();
	}

	// set radius
	void radiusIs(const double rad) {
		_graphic->setRadius(rad);
	}

	// get density
	double density() const {
		return _density;
	}

	// set density
	void densityIs(const double density) {
		if (density <= 0.0) {
			throw(std::runtime_error("Negative density."));
		}
		_density = density;
	}

	// get volume
	double volume() const {
		return pow(radius(),3)*4.0/3.0*M_PI;
	}

	// get mass
	double mass() const {
		return _density * volume();
	}

	// data members
public:
	// name of the fruit
	std::string _name;

	// graphics for this fruit
	chai3d::cShapeSphere* _graphic;

	// density of the fruit
	double _density;
};

#endif //FRUIT_H
