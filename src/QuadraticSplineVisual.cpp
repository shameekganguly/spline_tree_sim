#include "QuadraticSplineVisual.h"
#include <math.h>
#include <stdexcept>

#include <iostream>
using namespace std;

using namespace chai3d;
using namespace Eigen;

// ctor
QuadraticSplineVisual::QuadraticSplineVisual(const QuadraticSplineKinematic* spline_kinematic)
{
	_kinematic = spline_kinematic;

	// default values
	_nv_longitudinal = 12;
	_nv_plane = 20;

	// default material
    m_material = cMaterial::create();
    m_material->setShininess(10);
    m_material->m_ambient.set ((float)0.3, (float)0.3, (float)0.3);
    m_material->m_diffuse.set ((float)0.1, (float)0.7, (float)0.8);
    m_material->m_specular.set((float)1.0, (float)1.0, (float)1.0);
}

// dtor
QuadraticSplineVisual::~QuadraticSplineVisual() {
	// nothing to do
}

// render
// NOTE: this implementation is based on chai3d::cShapeBox
void QuadraticSplineVisual::render(cRenderOptions& a_options)
{
#ifdef C_USE_OPENGL

    /////////////////////////////////////////////////////////////////////////
    // ENABLE SHADER
    /////////////////////////////////////////////////////////////////////////
    if ((m_shaderProgram != nullptr) && (!a_options.m_creating_shadow_map))
    {
        // enable shader
        m_shaderProgram->use(this, a_options);
    }


    /////////////////////////////////////////////////////////////////////////
    // Render parts that use material properties
    /////////////////////////////////////////////////////////////////////////
    if (SECTION_RENDER_PARTS_WITH_MATERIALS(a_options, m_useTransparency))
    {
        // render material properties
        if (m_useMaterialProperty)
        {
            m_material->render(a_options);
        }


        if (!m_displayList.render(m_useDisplayList))
        {
            // create display list if requested
            m_displayList.begin(m_useDisplayList);

            // compute spacing parameters
            double _ds_longitudinal = _kinematic->_length/(_nv_longitudinal-1);
            double _ds_longitudinal_cap = _ds_longitudinal/3.0;
			double _ds_angle = 2.0*M_PI/_nv_plane;

            uint nv_endcap_offset = max((uint)10, uint(_kinematic->_radius/_ds_longitudinal_cap) + (uint)1);
            _ds_longitudinal_cap = _kinematic->_radius/(nv_endcap_offset - 1);
            uint nv_long_wendcaps = _nv_longitudinal + nv_endcap_offset*2;

            Vector3d point1, point2, point3, point4;
            //^ NOTE: point 1 and point 2 lie on the current plane
            // point 3 and 4 lie on the next plane
            // point 1 and 3 lie on the current angle
            // point 2 and 4 lie on the next angle
            Vector3d normal1, normal2;
            SplinePointPolar spoint;
            double s_i, s_ip1, t_i, t_ip1;
            // cout << "_ds_longitudinal " << _ds_longitudinal << " _ds_longitudinal_cap: " << _ds_longitudinal_cap << endl;
            // cout << "_nv_longitudinal " << _nv_longitudinal << " nv_endcap_offset: " << nv_endcap_offset << " nv_long_wendcaps: " << nv_long_wendcaps << endl;
            for (uint i = 0; i < (nv_long_wendcaps-1); i++) {
                // i = index of current cross section plane
                    // j = index of current vertex along the plane
                    if (i < nv_endcap_offset) {
                        s_i = ((double)i)*_ds_longitudinal_cap - _kinematic->_radius;
                    } else if (i > (nv_endcap_offset + _nv_longitudinal - 1)) {
                        s_i = ((double)(i - nv_endcap_offset - _nv_longitudinal))*_ds_longitudinal_cap + _kinematic->_length;
                    } else {
                        s_i = ((double)(i - nv_endcap_offset))*_ds_longitudinal;
                    }
                    if (i+1 < nv_endcap_offset) {
                        s_ip1 = ((double)i+1)*_ds_longitudinal_cap - _kinematic->_radius;
                    } else if (i+1 > (nv_endcap_offset + _nv_longitudinal - 1)) {
                        s_ip1 = ((double)(i+1 - nv_endcap_offset - _nv_longitudinal))*_ds_longitudinal_cap + _kinematic->_length;
                    } else {
                        s_ip1 = ((double)(i+1 - nv_endcap_offset))*_ds_longitudinal;
                    }
                    s_ip1 = min(s_ip1, _kinematic->_length + _kinematic->_radius);

                    // cout << "i: " << i << " s_i: " << s_i << " s_ip1 " << s_ip1 << endl;
                    t_i = _kinematic->radius(s_i);
                    t_ip1 = _kinematic->radius(s_ip1);

            	for (uint j = 0; j < (_nv_plane); j++) {
            		// get spatial locations of the four corners in the local
            		// frame. t != 0, eta = 0 => on the y-axis
                    spoint = SplinePointPolar(/* s */s_i, /* t */ t_i, /* eta */ ((double)j)*_ds_angle);
            		_kinematic->deformedLocation(point1, spoint);
                    spoint = SplinePointPolar(/* s */s_ip1, /* t */ t_ip1, /* eta */ ((double)j)*_ds_angle);
					_kinematic->deformedLocation(point3, spoint);
            		if (j == (_nv_plane - 1)) {
            			// wrap around. so point 2 and 4 are on zero angle again
                        spoint = SplinePointPolar(/* s */s_i, /* t */ t_i, /* eta */ 0.0);
						_kinematic->deformedLocation(point2, spoint);
                        spoint = SplinePointPolar(/* s */s_ip1, /* t */ t_ip1, /* eta */ 0.0);
						_kinematic->deformedLocation(point4, spoint);
            		} else {
                        spoint = SplinePointPolar(/* s */s_i, /* t */ t_i, /* eta */ ((double)j+1)*_ds_angle);
						_kinematic->deformedLocation(point2, spoint);
                        spoint = SplinePointPolar(/* s */s_ip1, /* t */ t_ip1, /* eta */ ((double)j+1)*_ds_angle);
						_kinematic->deformedLocation(point4, spoint);
            		}

            		// add triangle 1
            		// - compute normal
            		normal1 = (point3 - point1).cross(point2 - point1);
            		normal1.normalize(); // no fear of divide by zero here
            		// - add vertices
					glBegin(GL_POLYGON);
						glNormal3d(normal1[0], normal1[1], normal1[2]);
						glVertex3d(point1[0], point1[1], point1[2]);
						glVertex3d(point3[0], point3[1], point3[2]);
						glVertex3d(point2[0], point2[1], point2[2]);
					glEnd();

            		// add triangle 2
            		// - compute normal
            		normal2 = (point2 - point4).cross(point3 - point4);
            		normal2.normalize(); // no fear of divide by zero here
            		// - add vertices
            		glBegin(GL_POLYGON);
						glNormal3d(normal2[0], normal2[1], normal2[2]);
						glVertex3d(point2[0], point2[1], point2[2]);
						glVertex3d(point3[0], point3[1], point3[2]);
						glVertex3d(point4[0], point4[1], point4[2]);
					glEnd();
            	}
            }

            // finalize display list
            m_displayList.end(true);
        }
    }

    /////////////////////////////////////////////////////////////////////////
    // DISABLE SHADER
    /////////////////////////////////////////////////////////////////////////
    if ((m_shaderProgram != nullptr) && (!a_options.m_creating_shadow_map))
    {
        // disable shader
        m_shaderProgram->disable();
    }

#endif
}
