// 01-branch: visual display of a single branch in a tree

#include <Sai2Graphics.h>
#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew
#include "TreeKinematic.h"
#include "TreeVisual.h"
#include "timer/LoopTimer.h"

#include <Eigen/Core>
#include <iostream>
#include <string>
#include <thread>

using namespace std;
using namespace Eigen;
using namespace chai3d;

const string world_file = "resources/02-simple_tree/world.urdf";
const string camera_name = "camera_front";

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

// function for updating haptics
bool fSimulationRunning = false;
void updateHaptics(
	cGenericHapticDevicePtr device,
	chai3d::cShapeSphere* cursor,
	chai3d::cShapeLine* haptic_force_line,
	Sai2Graphics::Sai2Graphics* graphics,
	QuadraticSplineKinematic* spline,
	QuadraticSplineVisual* spline_graphic,
	chai3d::cShapeSphere* cherry,
	const double s_cherry
);

// flags for haptic interaction
bool fHapticDeviceEnabled = false;
bool fHapticSwitchPressed = false;

int main(int argc, char** argv) {
	cout << "Loading URDF world model file: " << world_file << endl;

	// load graphics scene
	auto graphics = new Sai2Graphics::Sai2Graphics(world_file, false);
	Vector3d camera_pos, camera_lookat, camera_vertical;
	graphics->getCameraPose(camera_name, camera_pos, camera_vertical, camera_lookat);
	graphics->_world->setBackgroundColor(0.3, 0.5, 0.7);

	// create a tree
	auto tree = new TreeKinematic("simple-tree");
	auto tree_visual = new TreeVisual(tree);
	Affine3d transform = tree->transformFromWorld();
	tree->transformFromWorldIs(transform.rotate(Quaterniond(0.707, 0.0, -0.707, 0.0)));
	auto trunk_br = tree->branch(tree->trunk());
	auto spline = trunk_br->spline();
	spline->_length = 0.9; spline->_radius = 0.07;

	// add a branch
	auto branch1 = tree->branchIs("branch-1", tree->trunk(), 0.0);
	ChildBranchInfo info;
	info.name = branch1->_name;
	info.s = 0.6;
	info.rotation << 0.866, -0.5, 0.0,
					0.5, 0.866, 0.0,
					0.0, 0.0, 1.0;
	tree->branchInfoInParentUpdate(branch1->_name, info);
	spline = branch1->spline();
	spline->_length = 0.8; spline->_radius = 0.05;

	// add a branch
	auto branch2 = tree->branchIs("branch-2", tree->trunk(), 0.0);
	info.name = branch2->_name;
	info.s = 0.8;
	info.rotation << 0.866, 0.5, 0.0,
					-0.5, 0.866, 0.0,
					0.0, 0.0, 1.0;
	tree->branchInfoInParentUpdate(branch2->_name, info);
	spline = branch2->spline();
	spline->_length = 0.6; spline->_radius = 0.04;

	// add a branch
	auto branch3 = tree->branchIs("branch-3", tree->trunk(), 0.9);
	spline = branch3->spline();
	spline->_length = 0.3; spline->_radius = 0.02;

	// create the cherry
	double r0 = branch1->spline()->_radius * 2;
	double s_cherry = branch1->spline()->_length;
	auto cherry_fruit = tree->fruitIs("cherry", branch1->_name, s_cherry);
	cherry_fruit->radiusIs(r0);
	auto cherry = cherry_fruit->graphic();
	cherry->m_material->setBrownMaroon();
	cherry->m_material->setShininess(100);

	// inspect cherry position
	Vector3d cherry_pos;
	tree->positionInWorld(cherry_pos, branch1->_name, s_cherry);
	cout << "Cherry located at: " << cherry_pos.transpose() << endl;

	// change material for branches and fruits
	auto material = trunk_br->splineVisual()->m_material;
	material->setColorf(0.3, 0.15, 0.1);
	material->setShininess(100);
	tree_visual->branchMaterialIs(material);

	// inspect Jacobian to cherry
	MatrixXd cherry_Jv;
	tree->jacobianLinear(cherry_Jv, branch1->_name, s_cherry);
	cout << "Cherry Jacobian: " << endl;
	cout << cherry_Jv << endl;
	cout << "Branch 0: " << tree->branchAtIndex(0) << endl;
	cout << "Branch 1: " << tree->branchAtIndex(1) << endl;
	cout << "Branch 2: " << tree->branchAtIndex(2) << endl;
	cout << "Branch 3: " << tree->branchAtIndex(3) << endl;

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

	/*------- Set up haptic interaction -------*/
	// create a sphere (cursor) to represent the haptic device
    auto cursor = new cShapeSphere(0.1);
    // create a line to show the haptic force
    auto haptic_force_line = new cShapeLine();
    haptic_force_line->setShowEnabled(false);
    haptic_force_line->setLineWidth(4.0);
		
	// create a haptic device handler
    auto handler = new cHapticDeviceHandler();

	// get a handle to the first haptic device
    cGenericHapticDevicePtr hapticDevice;
	if (!handler->getDevice(hapticDevice, 0)) {
		cout << "No haptic device found. " << endl;
		fHapticDeviceEnabled = false;
	} else {
		hapticDevice->open();
		hapticDevice->calibrate();
		fHapticDeviceEnabled = true;
		// if the device has a gripper, enable the gripper to simulate a user switch
	    hapticDevice->setEnableGripperUserSwitch(true);
	    graphics->_world->addChild(cursor);
	    graphics->_world->addChild(haptic_force_line);
	}

	// thread haptics_thread(updateHaptics, hapticDevice, cursor, haptic_force_line, graphics, spline, spline_graphic, cherry, s_cherry);

	/*------- Loop -------*/
	// cache variables
	double last_cursorx, last_cursory;

	Eigen::MatrixXd G;
	Eigen::Matrix3d R;
	Eigen::Vector3d center_point = Eigen::Vector3d::Zero();

	// while window is open:
    while (!glfwWindowShouldClose(window))
	{
		// update graphics. this automatically waits for the correct amount of time
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		tree_visual->updateGraphics(graphics->_world);

		// render scene
		graphics->_world->updateShadowMaps(false);
		graphics->render(camera_name, width, height);
		
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

	// stop simulation
	fSimulationRunning = false;
	// haptics_thread.join();
	
    // destroy context
    glfwDestroyWindow(window);

    // terminate
    glfwTerminate();

	return 0;
}

//------------------------------------------------------------------------------
void updateHaptics(
	cGenericHapticDevicePtr hapticDevice,
	chai3d::cShapeSphere* cursor,
	chai3d::cShapeLine* haptic_force_line,
	Sai2Graphics::Sai2Graphics* graphics,
	QuadraticSplineKinematic* spline,
	QuadraticSplineVisual* spline_graphic,
	chai3d::cShapeSphere* cherry,
	const double s_cherry
) {
	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(1000); //1000Hz timer
	double last_time = timer.elapsedTime(); //secs

	bool fTimerDidSleep = true;

	// cherry updates
	Vector3d cherry_pos_local, cherry_pos_global;
	bool cherry_still_on_plant = true;
	const double cherry_pluck_force_thresh = 1.5;
	Vector3d cherry_vel, cherry_acc, cherry_force;

	// spline dynamics variables
	double ks = 10.0;
	double b = 0.05;
	double cherry_r = cherry->getRadius();
	double cherry_r_max = 0.15;
	const double cherry_growth_rate = 0.007; // r/ sec
	MatrixXd Jv_s;
	MatrixXd Jlp;
	VectorXd dq(2);
	Vector3d F_cherry;
	VectorXd gamma_cherry;
	const double density_cherry = 0.65;

	// haptics device
	cVector3d position;
	const double scale_factor = 80.0;
	Vector3d home_pos(-0.25, 0.0, 1.5);
	Vector3d F_haptic;
	double last_distance;
	Vector3d last_dir;
	const double haptic_kp = 1.5;
	const double haptic_force_scale = 9.0;

	// start simulation loop
	fSimulationRunning = true;
	while (fSimulationRunning) { //automatically set to false when simulation is quit
		fTimerDidSleep = timer.waitForNextLoop();

		// update time
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time;

		// /* --- SPLINE KINEMATICS UPDATES BEG --- */
		// if (cherry_still_on_plant) {
		// 	spline->splineLocation(cherry_pos_local, s_cherry);
		// 	cherry->setLocalPos(cherry_pos_local);

		// 	if (cherry_r < cherry_r_max) {
		// 		cherry_r += cherry_growth_rate*loop_dt;
		// 	}
		// 	cherry->setRadius(cherry_r);
		// }
		// /* --- SPLINE KINEMATICS UPDATES END --- */

		// if (fHapticDeviceEnabled) {
		// 	// haptics updates
	 //        hapticDevice->getPosition(position);
	 //        cursor->setLocalPos(position*scale_factor + home_pos);

	 //        if (cherry_still_on_plant || fHapticSwitchPressed) {
		//         hapticDevice->getUserSwitch(0, fHapticSwitchPressed);
	 //        }

	 //        if (fHapticSwitchPressed) {
	 //        	// ignore if the cursor is too far
	 //        	if (
	 //        		(cursor->getLocalPos().eigen() - cherry->getGlobalPos().eigen()).norm() > (cherry->getRadius()*1.0 + cursor->getRadius())
	 //        		&& !haptic_force_line->getShowEnabled()
  //       		) {
	 //        		fHapticSwitchPressed = false;
	 //        	} else {
		//         	haptic_force_line->m_pointA = cherry->getGlobalPos();
		//         	haptic_force_line->m_pointB = cursor->getLocalPos();
		//         	haptic_force_line->setShowEnabled(true);

		//         	// compute forces
		//         	if (cherry_still_on_plant) {
		//         		F_haptic = haptic_kp * (haptic_force_line->m_pointB - haptic_force_line->m_pointA).eigen();
		//         	} else {
		//         		F_haptic.setZero();
		//         		F_haptic[2] = 9.8*density_cherry*4.0/3.0*M_PI*pow(cherry_r,3);
		//         	}
	 //        	}
		// 	}
		// 	if (!fHapticSwitchPressed) {
		// 		F_haptic.setZero();
		// 		haptic_force_line->setShowEnabled(false);
	 //        }
		// }

		// /* --- CHERRY DYNAMICS BEG ---*/
		// if (F_haptic.norm() > cherry_pluck_force_thresh) {
		// 	// remove cherry from plant and add it to tree
		// 	cherry_pos_global = cherry->getGlobalPos().eigen();
		// 	spline_graphic->removeChild(cherry);
		// 	graphics->_world->addChild(cherry);
		// 	// set flag
		// 	cherry_still_on_plant = false;
		// }

		// if (!cherry_still_on_plant) {
		// 	if (fHapticSwitchPressed) {
		// 		cherry->setLocalPos(cursor->getLocalPos());
		// 		cherry_vel.setZero();
		// 	} else {
		// 		cherry_pos_global = cherry->getGlobalPos().eigen();
		// 		cherry_acc.setZero();
		// 		cherry_acc[2] -= 9.8;
		// 		cherry_acc += 1.0/(density_cherry*4.0/3.0*M_PI*pow(cherry_r,3))*F_haptic;
		// 		cherry_vel += cherry_acc*0.01;
		// 		cherry_pos_global += cherry_vel*0.01;
		// 		// clamp to ground
		// 		double cherry_min_height_ground = 0.0;
		// 		if (cherry_pos_global[2] < cherry_min_height_ground) {
		// 			cherry_pos_global[2] = cherry_min_height_ground;
		// 		}
		// 		cherry->setLocalPos(cherry_pos_global);
		// 	}
		// }

		// /* --- CHERRY DYNAMICS END ---*/

		// /* --- SPLINE DYNAMICS BEG ---*/
		// spline->splineLinearJacobian(Jv_s, s_cherry);
		// F_cherry.setZero();
		// if (cherry_still_on_plant) {
		// 	F_cherry[2] = -9.8*density_cherry*4.0/3.0*M_PI*pow(cherry_r,3);
		// 	F_cherry += F_haptic;
		// }
		// gamma_cherry = Jv_s.transpose()*F_cherry;

		// spline->splineProjectionLengthJacobian(Jlp, 0.5);
		// dq[0] = -Jlp(0,0)*ks/b*spline->splineProjectionLength(0.5) + gamma_cherry[0]/b;
		// dq[1] = -Jlp(0,1)*ks/b*spline->splineProjectionLength(0.5) + gamma_cherry[1]/b;
		// spline->_alpha += dq[0]*loop_dt;
		// spline->_beta += dq[1]*loop_dt;

		// if (fHapticDeviceEnabled) {
		// 	// apply haptics forces
		// 	hapticDevice->setForceAndTorqueAndGripperForce(-cVector3d(F_haptic)*haptic_force_scale, Vector3d(0.0, 0.0, 0.0), 0.0);
		// }

		// -------------------------------------------
		// update last time
		last_time = curr_time;
	}
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
