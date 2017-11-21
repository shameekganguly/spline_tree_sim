// 01-branch: visual display of a single branch in a tree

#include <Sai2Graphics.h>
#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew
#include "QuadraticSplineKinematic.h"
#include "QuadraticSplineVisual.h"

#include <Eigen/Core>
#include <iostream>
#include <string>

using namespace std;
using namespace Eigen;

const string world_file = "resources/01-branch/world.urdf";
const string camera_name = "camera_isometric";

// callback to print glfw errors
void glfwError(int error, const char* description);

// callback when a key is pressed
void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods);

// callback when a mouse button is pressed
void mouseClick(GLFWwindow* window, int button, int action, int mods);

// flags for scene camera movement
bool fTransXp = false;
bool fTransXn = false;
bool fTransYp = false;
bool fTransYn = false;
bool fRotPanTilt = false;

int main(int argc, char** argv) {
	cout << "Loading URDF world model file: " << world_file << endl;

	// load graphics scene
	auto graphics = new Sai2Graphics::Sai2Graphics(world_file, false);
	Vector3d camera_pos, camera_lookat, camera_vertical;
	graphics->getCameraPose(camera_name, camera_pos, camera_vertical, camera_lookat);

	// create the kinematic spline element
	// TODO: store the transform from the world frame to the local frame of
	// this spline
	auto spline = new QuadraticSplineKinematic();
	spline->_length = 1.0;
	spline->_radius = 0.04;
	spline->_alpha = 0.00001;
	spline->_beta = 0.00001;

	// create the spline element to display
	auto spline_graphic = new QuadraticSplineVisual(spline);
	graphics->_world->addChild(spline_graphic);
	spline_graphic->setLocalPos(Vector3d(-0.25, 0.0, 1.5));
	spline_graphic->m_material->setColorf(0.3, 0.15, 0.1);
	spline_graphic->m_material->setShininess(100);
	spline_graphic->_nv_longitudinal = 50;

	// create the cherry
	double r0 = spline->_radius;
	double s_cherry = spline->_length;
	auto cherry = new chai3d::cShapeSphere(r0);
	spline_graphic->addChild(cherry);
	cherry->m_material->setBrownMaroon();
	cherry->m_material->setShininess(100);

	/*------- Set up visualization -------*/
    // set up error callback
    glfwSetErrorCallback(glfwError);

    // initialize GLFW
    glfwInit();

    // retrieve resolution of computer display and position window accordingly
    GLFWmonitor* primary = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode = glfwGetVideoMode(primary);

    // information about computer screen and GLUT display window
	int screenW = mode->width;
    int screenH = mode->height;
    int windowW = 0.8 * screenH;
    int windowH = 0.5 * screenH;
    int windowPosY = (screenH - windowH) / 2;
    int windowPosX = windowPosY;

    // create window and make it current
    glfwWindowHint(GLFW_VISIBLE, 0);
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "01-branch", NULL, NULL);
	glfwSetWindowPos(window, windowPosX, windowPosY);
	glfwShowWindow(window);
    glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

    // set callbacks
	glfwSetKeyCallback(window, keySelect);
	glfwSetMouseButtonCallback(window, mouseClick);

	// cache variables
	double last_cursorx, last_cursory;

	Eigen::MatrixXd G;
	Eigen::Matrix3d R;
	Eigen::Vector3d center_point = Eigen::Vector3d::Zero();

	// cherry position
	Vector3d cherry_pos_local;

	// spline dynamics variables
	double ks = 0.15;
	double b = 0.05;
	double dt = 0.01;
	double cherry_r = r0;
	double cherry_r_max = 0.25;
	const double cherry_growth_rate = 0.02; // r/ sec
	MatrixXd Jv_s;
	MatrixXd Jlp;
	VectorXd dq(2);
	Vector3d F_cherry;
	VectorXd gamma_cherry;
	const double density_cherry = 0.15;

	// while window is open:
    while (!glfwWindowShouldClose(window))
	{
		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		
		/* --- SPLINE KINEMATICS UPDATES BEG --- */
		spline->splineLocation(cherry_pos_local, s_cherry);
		cherry->setLocalPos(cherry_pos_local);

		if (cherry_r < cherry_r_max) {
			cherry_r += cherry_growth_rate*dt;
		}
		cherry->setRadius(cherry_r);
		/* --- SPLINE KINEMATICS UPDATES END --- */

		graphics->render(camera_name, width, height);


		/* --- SPLINE DYNAMICS BEG ---*/
		spline->splineLinearJacobian(Jv_s, s_cherry);
		F_cherry.setZero();
		F_cherry[2] = -9.8*density_cherry*4.0/3.0*M_PI*pow(cherry_r,3);
		gamma_cherry = Jv_s.transpose()*F_cherry;

		spline->splineTipProjectionLengthJacobian(Jlp);
		dq[0] = -Jlp(0,0)*ks/b*spline->splineTipProjectionLength() + gamma_cherry[0]/b;
		dq[1] = -Jlp(0,1)*ks/b*spline->splineTipProjectionLength() + gamma_cherry[1]/b;
		spline->_alpha += dq[0]*dt;
		spline->_beta += dq[1]*dt;
		
		/* --- SPLINE DYNAMICS END ---*/

		// swap buffers
		glfwSwapBuffers(window);

		// wait until all GL commands are completed
		glFinish();

		// check for any OpenGL errors
		GLenum err;
		err = glGetError();
		assert(err == GL_NO_ERROR);

	    // poll for events
	    glfwPollEvents();
	
		// move scene camera as required
    	// graphics->getCameraPose(camera_name, camera_pos, camera_vertical, camera_lookat);
    	Eigen::Vector3d cam_up_axis;
    	// cam_up_axis = camera_vertical;
    	// cam_up_axis.normalize();
    	cam_up_axis << 0.0, 0.0, 1.0; //TODO: there might be a better way to do this
	    Eigen::Vector3d cam_roll_axis = (camera_lookat - camera_pos).cross(cam_up_axis);
    	cam_roll_axis.normalize();
    	Eigen::Vector3d cam_lookat_axis = camera_lookat;
    	cam_lookat_axis.normalize();
    	if (fTransXp) {
	    	camera_pos = camera_pos + 0.05*cam_roll_axis;
	    	camera_lookat = camera_lookat + 0.05*cam_roll_axis;
	    }
	    if (fTransXn) {
	    	camera_pos = camera_pos - 0.05*cam_roll_axis;
	    	camera_lookat = camera_lookat - 0.05*cam_roll_axis;
	    }
	    if (fTransYp) {
	    	// camera_pos = camera_pos + 0.05*cam_lookat_axis;
	    	camera_pos = camera_pos + 0.05*cam_up_axis;
	    	camera_lookat = camera_lookat + 0.05*cam_up_axis;
	    }
	    if (fTransYn) {
	    	// camera_pos = camera_pos - 0.05*cam_lookat_axis;
	    	camera_pos = camera_pos - 0.05*cam_up_axis;
	    	camera_lookat = camera_lookat - 0.05*cam_up_axis;
	    }
	    if (fRotPanTilt) {
	    	// get current cursor position
	    	double cursorx, cursory;
			glfwGetCursorPos(window, &cursorx, &cursory);
			//TODO: might need to re-scale from screen units to physical units
			double compass = 0.006*(cursorx - last_cursorx);
			double azimuth = 0.006*(cursory - last_cursory);
			double radius = (camera_pos - camera_lookat).norm();
			Eigen::Matrix3d m_tilt; m_tilt = Eigen::AngleAxisd(azimuth, -cam_roll_axis);
			camera_pos = camera_lookat + m_tilt*(camera_pos - camera_lookat);
			Eigen::Matrix3d m_pan; m_pan = Eigen::AngleAxisd(compass, -cam_up_axis);
			camera_pos = camera_lookat + m_pan*(camera_pos - camera_lookat);
	    }
	    graphics->setCameraPose(camera_name, camera_pos, cam_up_axis, camera_lookat);
	    glfwGetCursorPos(window, &last_cursorx, &last_cursory);

	}

    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

//------------------------------------------------------------------------------

void glfwError(int error, const char* description) {
	cerr << "GLFW Error: " << description << endl;
	exit(1);
}

//------------------------------------------------------------------------------

void keySelect(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	bool set = (action != GLFW_RELEASE);
    switch(key) {
		case GLFW_KEY_ESCAPE:
			// exit application
			glfwSetWindowShouldClose(window,GL_TRUE);
			break;
		case GLFW_KEY_RIGHT:
			fTransXp = set;
			break;
		case GLFW_KEY_LEFT:
			fTransXn = set;
			break;
		case GLFW_KEY_UP:
			fTransYp = set;
			break;
		case GLFW_KEY_DOWN:
			fTransYn = set;
			break;
		default:
			break;
    }
}

//------------------------------------------------------------------------------

void mouseClick(GLFWwindow* window, int button, int action, int mods) {
	bool set = (action != GLFW_RELEASE);
	//TODO: mouse interaction with robot
	switch (button) {
		// left click pans and tilts
		case GLFW_MOUSE_BUTTON_LEFT:
			fRotPanTilt = set;
			// NOTE: the code below is recommended but doesn't work well
			// if (fRotPanTilt) {
			// 	// lock cursor
			// 	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
			// } else {
			// 	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			// }
			break;
		// if right click: don't handle. this is for menu selection
		case GLFW_MOUSE_BUTTON_RIGHT:
			//TODO: menu
			break;
		// if middle click: don't handle. doesn't work well on laptops
		case GLFW_MOUSE_BUTTON_MIDDLE:
			break;
		default:
			break;
	}
}
