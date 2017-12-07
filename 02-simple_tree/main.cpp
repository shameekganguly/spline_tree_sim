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
	TreeKinematic* tree
);
string getClosestFruitToPoint(const TreeKinematic* tree, const Vector3d point_pos);

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

	// create a cherry
	double r0 = branch1->spline()->_radius * 2;
	double s_cherry = branch1->spline()->_length;
	auto cherry = tree->fruitIs("cherry", branch1->_name, s_cherry);
	cherry->radiusIs(r0);
	cherry->densityIs(10);

	// inspect cherry position
	Vector3d cherry_pos;
	tree->positionInWorld(cherry_pos, branch1->_name, s_cherry);
	cout << "Cherry located at: " << cherry_pos.transpose() << endl;

	// inspect Jacobian to cherry
	MatrixXd cherry_Jv;
	tree->jacobianLinear(cherry_Jv, branch1->_name, s_cherry);
	cout << "Cherry Jacobian: " << endl;
	cout << cherry_Jv << endl;
	cout << "Branch 0: " << tree->branchAtIndex(0) << endl;
	cout << "Branch 1: " << tree->branchAtIndex(1) << endl;
	cout << "Branch 2: " << tree->branchAtIndex(2) << endl;
	cout << "Branch 3: " << tree->branchAtIndex(3) << endl;

	// create another cherry
	double r02 = branch2->spline()->_radius * 2;
	double s_cherry2 = branch2->spline()->_length;
	auto cherry2 = tree->fruitIs("cherry2", branch2->_name, s_cherry2);
	cherry2->radiusIs(r02);
	cherry2->densityIs(10);

	// change material for branches and fruits
	auto material = trunk_br->splineVisual()->m_material;
	material->m_diffuse = cColorf(0.3, 0.15, 0.1);
	material->m_ambient = cColorf(0.1, 0.06, 0.0);
	material->m_specular = cColorf(0.0, 0.05, 0.05);
	material->setShininess(100);
	tree_visual->branchMaterialIs(material);

	auto fruit_material = cherry->graphic()->m_material;
	fruit_material->m_diffuse = cColorf(0.6, 0.4, 0.05);
	fruit_material->m_ambient = cColorf(0.2, 0.02, 0.02);
	fruit_material->m_specular = cColorf(0.0, 0.05, 0.05);
	fruit_material->setShininess(100);
	tree_visual->fruitMaterialIs(fruit_material);


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
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "02-simple_tree", NULL, NULL);
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

	thread haptics_thread(updateHaptics, hapticDevice, cursor, haptic_force_line, graphics, tree);

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
	haptics_thread.join();
	
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
	TreeKinematic* tree
) {
	// create a timer
	LoopTimer timer;
	timer.initializeTimer();
	timer.setLoopFrequency(1000); //1000Hz timer
	double last_time = timer.elapsedTime(); //secs

	bool fTimerDidSleep = true;

	// cherry updates
	Fruit* cherry;
	Vector3d cherry_pos_local, cherry_pos_global;
	bool cherry_still_on_plant = true;
	const double cherry_pluck_force_thresh = 1.5;
	Vector3d cherry_acc;

	// spline dynamics variables
	double ks = 10000.0;
	double b = 0.05;
	MatrixXd Jv_cherry, Jv_haptic;
	MatrixXd Jlp;
	VectorXd dq(2);
	Vector3d F_cherry;
	VectorXd gamma_cherry(tree->dof());
	string cherry_branch_name;
	FruitInfo cherry_info;
	TreeKinematic::FruitList::iterator fruit_itr;
	TreeKinematic::BranchList::iterator branch_itr;
	BranchKinematic* branch_ptr;
	QuadraticSplineKinematic* spline_ptr;
	TreeKinematic::FruitList fallen_cherries;
	map<string, Vector3d> falling_cherry_vels;

	// haptics device
	cVector3d position;
	const double scale_factor = 80.0;
	Vector3d home_pos(-0.25, 0.0, 1.5);
	Vector3d F_haptic;
	double cursor_distance;
	Vector3d last_dir;
	const double haptic_kp = 1.5;
	const double haptic_force_scale = 9.0;

	bool carrying_cherry = false;
	string closest_cherry_name;
	string closest_cherry_branch_name;
	FruitInfo closest_cherry_info;
	Fruit* haptic_cherry;
	Vector3d closest_cherry_pos;

	// start simulation loop
	fSimulationRunning = true;
	while (fSimulationRunning) { //automatically set to false when simulation is quit
		fTimerDidSleep = timer.waitForNextLoop();

		// update time
		double curr_time = timer.elapsedTime();
		double loop_dt = curr_time - last_time;

		if (fHapticDeviceEnabled) {
			// haptics updates
	        hapticDevice->getPosition(position);
	        cursor->setLocalPos(position*scale_factor + home_pos);

	        hapticDevice->getUserSwitch(0, fHapticSwitchPressed);

	        if (fHapticSwitchPressed) {
	        	if (!carrying_cherry) {
	        		cherry_still_on_plant = true;
					// update cursor distance to closest fruit on tree
					closest_cherry_name = getClosestFruitToPoint(tree, cursor->getLocalPos().eigen());
					haptic_cherry = tree->fruit(closest_cherry_name);
					if (NULL == haptic_cherry) {
						fHapticSwitchPressed = false;
					} else {
						closest_cherry_branch_name = tree->fruitParentBranch(closest_cherry_name);
						closest_cherry_info = tree->branch(closest_cherry_branch_name)->fruitInfo(closest_cherry_name);
						tree->positionInWorld(closest_cherry_pos, closest_cherry_branch_name, closest_cherry_info.s);
						cursor_distance = (cursor->getLocalPos().eigen() - closest_cherry_pos).norm();
						// cout << "Cursor distance " << cursor_distance << endl;
						if (cursor_distance > (haptic_cherry->radius()*1.0 + cursor->getRadius())) {
							fHapticSwitchPressed = false;
						}
					}
	        	}
	        }

	        if(fHapticSwitchPressed) {
	        	carrying_cherry = true;
				haptic_force_line->m_pointA = haptic_cherry->graphic()->getLocalPos();
				haptic_force_line->m_pointB = cursor->getLocalPos();
				haptic_force_line->setShowEnabled(true);

				// compute forces
				if (cherry_still_on_plant) {
					F_haptic = haptic_kp * (haptic_force_line->m_pointB - haptic_force_line->m_pointA).eigen();
				} else {
					F_haptic.setZero();
					F_haptic[2] = 9.8*haptic_cherry->mass();
				}
	        }

			if (!fHapticSwitchPressed) {
				F_haptic.setZero();
				haptic_force_line->setShowEnabled(false);
				carrying_cherry = false;
	        }
		} else {
			F_haptic.setZero();
		}

		/* --- CHERRY DYNAMICS BEG ---*/
		if (F_haptic.norm() > cherry_pluck_force_thresh) {
			// remove cherry from plant and add it to tree
			tree->positionInWorld(closest_cherry_pos, closest_cherry_branch_name, closest_cherry_info.s);
			tree->fruitRem(haptic_cherry->_name);

			// TODO: technically this is not necessary since fruitRem does not
			// remove the visual element from the graphics scene
			graphics->_world->addChild(haptic_cherry->graphic());

			// add cherry to list of fallen fruits
			fallen_cherries[haptic_cherry->_name] = haptic_cherry;
			falling_cherry_vels[haptic_cherry->_name].setZero();
			cherry_still_on_plant = false;
		}

		// loop over fallen cherries
		for (auto fallen_fruit_itr: fallen_cherries) {
			// is fruit the same as the current haptic fruit
			if (fHapticSwitchPressed && fallen_fruit_itr.first.compare(haptic_cherry->_name) == 0) {
				fallen_fruit_itr.second->graphic()->setLocalPos(cursor->getLocalPos());
				falling_cherry_vels[fallen_fruit_itr.first].setZero();
			} else {
				cherry = fallen_fruit_itr.second;
				cherry_pos_global = cherry->graphic()->getLocalPos().eigen();
				cherry_acc.setZero();
				cherry_acc[2] -= 9.8;
				falling_cherry_vels[fallen_fruit_itr.first] += cherry_acc*0.01;
				cherry_pos_global += falling_cherry_vels[fallen_fruit_itr.first]*0.01;
				// clamp to ground
				double cherry_min_height_ground = 0.0;
				if (cherry_pos_global[2] < cherry_min_height_ground) {
					cherry_pos_global[2] = cherry_min_height_ground;
				}
				cherry->graphic()->setLocalPos(cherry_pos_global);
			}
		}
		/* --- CHERRY DYNAMICS END ---*/

		/* --- TREE DYNAMICS BEG ---*/
		gamma_cherry.setZero();
		// loop over fruits
		for (fruit_itr=tree->fruitsItrBegin(); fruit_itr!=tree->fruitsItrEnd(); ++fruit_itr) {
			cherry = fruit_itr->second;
			cherry_branch_name = tree->fruitParentBranch(fruit_itr->first);
			double s_cherry = tree->branch(cherry_branch_name)->fruitInfo(fruit_itr->first).s;

			F_cherry.setZero();
			F_cherry[2] = -9.8*fruit_itr->second->mass();
			if (carrying_cherry && fruit_itr->first.compare(haptic_cherry->_name) == 0) {
				F_cherry += F_haptic;
			}

			tree->jacobianLinear(Jv_cherry, cherry_branch_name, s_cherry);
			gamma_cherry += Jv_cherry.transpose()*F_cherry;
		}

		// cout << "Gamma cherry" << gamma_cherry.transpose() << endl;
		// loop over all branches
		for (branch_itr=tree->branchesItrBegin(); branch_itr!=tree->branchesItrEnd(); ++branch_itr) {
			branch_ptr = branch_itr->second;
			spline_ptr = branch_ptr->spline();
			// get index in gamma
			uint branch_index = tree->branchIndex(branch_ptr->_name);
			// compute dq
			spline_ptr->splineProjectionLengthJacobian(Jlp, 0.5);
			double ks_spline = ks*pow(spline_ptr->_radius, 2);
			dq[0] = -Jlp(0,0)*ks_spline/b*spline_ptr->splineProjectionLength(0.5) + gamma_cherry[2*branch_index + 0]/b;
			dq[1] = -Jlp(0,1)*ks_spline/b*spline_ptr->splineProjectionLength(0.5) + gamma_cherry[2*branch_index + 1]/b;
			spline_ptr->_alpha += dq[0]*loop_dt;
			spline_ptr->_beta += dq[1]*loop_dt;
		}
		/* --- TREE DYNAMICS END ---*/

		if (fHapticDeviceEnabled) {
			// apply haptics forces
			hapticDevice->setForceAndTorqueAndGripperForce(-cVector3d(F_haptic)*haptic_force_scale, Vector3d(0.0, 0.0, 0.0), 0.0);
		}

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

//------------------------------------------------------------------------------
string getClosestFruitToPoint(const TreeKinematic* tree, const Vector3d point_pos) {
	double min_distance = -1;
	string closest_fruit = "";

	// temp variables
	string branch_name;
	string fruit_name;
	const BranchKinematic* branch_ptr;
	FruitInfo fruit_info;
	Vector3d fruit_pos;
	double distance;
	for (auto fruit_itr = tree->fruitsItrBegin(); fruit_itr != tree->fruitsItrEnd(); ++fruit_itr) {
		fruit_name = fruit_itr->first;
		branch_name = tree->fruitParentBranch(fruit_name);
		branch_ptr = tree->branch(branch_name);
		fruit_info = branch_ptr->fruitInfo(fruit_name);
		tree->positionInWorld(fruit_pos, branch_name, fruit_info.s);
		distance = (fruit_pos - point_pos).norm();
		if (min_distance < 0 || min_distance > distance) {
			min_distance = distance;
			closest_fruit = fruit_name;
		}
	}
	return closest_fruit;
}
