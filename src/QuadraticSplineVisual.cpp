#include "QuadraticSplineVisual.h"
#include <math.h>

using namespace chai3d;

// ctor
QuadraticSplineVisual::QuadraticSplineVisual()
{
	// default values
	_radius = 0.03;
	_length = 0.1;
	_alpha = 0.0;
	_beta = 0.0;
	_nv_longitudinal = 4;
	_ds_longitudinal = _length/_nv_longitudinal;
	_nv_plane = 6;
	_ds_angle = M_PI/_nv_plane;
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
            // // create display list if requested
            // m_displayList.begin(m_useDisplayList);

            // // render box
            // glBegin(GL_POLYGON);
            //     glNormal3d(1.0, 0.0, 0.0);
            //     glVertex3d(m_hSizeX, m_hSizeY, m_hSizeZ);
            //     glVertex3d(m_hSizeX,-m_hSizeY, m_hSizeZ);
            //     glVertex3d(m_hSizeX,-m_hSizeY,-m_hSizeZ);
            //     glVertex3d(m_hSizeX, m_hSizeY,-m_hSizeZ);
            // glEnd();
        
            // glBegin(GL_POLYGON);
            //     glNormal3d(-1.0, 0.0, 0.0);
            //     glVertex3d(-m_hSizeX, m_hSizeY,-m_hSizeZ);
            //     glVertex3d(-m_hSizeX,-m_hSizeY,-m_hSizeZ);
            //     glVertex3d(-m_hSizeX,-m_hSizeY, m_hSizeZ);
            //     glVertex3d(-m_hSizeX, m_hSizeY, m_hSizeZ);
            // glEnd();

            // glBegin(GL_POLYGON);
            //     glNormal3d(0.0, 1.0, 0.0);
            //     glVertex3d( m_hSizeX, m_hSizeY,-m_hSizeZ);
            //     glVertex3d(-m_hSizeX, m_hSizeY,-m_hSizeZ);
            //     glVertex3d(-m_hSizeX, m_hSizeY, m_hSizeZ);
            //     glVertex3d( m_hSizeX, m_hSizeY, m_hSizeZ);
            // glEnd();
        
            // glBegin(GL_POLYGON);
            //     glNormal3d(0.0,-1.0, 0.0);
            //     glVertex3d( m_hSizeX,-m_hSizeY, m_hSizeZ);
            //     glVertex3d(-m_hSizeX,-m_hSizeY, m_hSizeZ);
            //     glVertex3d(-m_hSizeX,-m_hSizeY,-m_hSizeZ);
            //     glVertex3d( m_hSizeX,-m_hSizeY,-m_hSizeZ);
            // glEnd();

            // glBegin(GL_POLYGON);
            //     glNormal3d(0.0, 0.0, 1.0);
            //     glVertex3d( m_hSizeX, m_hSizeY, m_hSizeZ);
            //     glVertex3d(-m_hSizeX, m_hSizeY, m_hSizeZ);
            //     glVertex3d(-m_hSizeX,-m_hSizeY, m_hSizeZ);
            //     glVertex3d( m_hSizeX,-m_hSizeY, m_hSizeZ);
            // glEnd();
            // glBegin(GL_POLYGON);
            //     glNormal3d(0.0, 0.0,-1.0);
            //     glVertex3d( m_hSizeX,-m_hSizeY,-m_hSizeZ);
            //     glVertex3d(-m_hSizeX,-m_hSizeY,-m_hSizeZ);
            //     glVertex3d(-m_hSizeX, m_hSizeY,-m_hSizeZ);
            //     glVertex3d( m_hSizeX, m_hSizeY,-m_hSizeZ);
            // glEnd();

            // // finalize display list
            // m_displayList.end(true);
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
