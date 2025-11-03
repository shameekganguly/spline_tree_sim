#include <iostream>
#include <string>
#include <string_view>
#include <thread>

// clang-format off
#include <chai3d.h>
#include <GLFW/glfw3.h> //must be loaded after loading opengl/glew
// clang-format on
#include <Eigen/Core>

#include "QuadraticSplineKinematic.h"
#include "QuadraticSplineVisual.h"
#include "TreeDynamics.h"
#include "TreeKinematic.h"
#include "TreeParser.h"
#include "TreeVisual.h"
#include "graphics/Graphics.h"
#include "timer/LoopTimer.h"

using namespace Eigen;
using namespace chai3d;

// callback to print glfw errors
void glfwError(int error, const char *description);

// callback when a key is pressed
void keySelect(GLFWwindow *window, int key, int scancode, int action, int mods);

// callback when a mouse button is pressed
void mouseClick(GLFWwindow *window, int button, int action, int mods);

// flags for scene camera movement
bool fTransXp = false;
bool fTransXn = false;
bool fTransYp = false;
bool fTransYn = false;
bool fRotPanTilt = false;

// function for updating scene
bool fSimulationRunning = false;
void update(TreeKinematic *tree_kinematic);

const std::string kCameraName = "camera";

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Did not pass model file" << std::endl;
    return 0;
  }
  std::string model_file(argv[1]);
  std::cout << "Loading: " << model_file << std::endl;
  auto graphics = std::make_unique<spline_sim::Graphics>();

  // Add lights
  graphics->CreateLight(/*pos=*/cVector3d(3, -3, 3),
                        /*look_at=*/cVector3d(0, 0, 2));
  graphics->CreateLight(/*pos=*/cVector3d(3, 3, 3),
                        /*look_at=*/cVector3d(0, 0, 2));

  // initialize a chai camera
  cCamera *camera = graphics->CreateCamera(kCameraName);

  // position and orient the camera
  Vector3d camera_pos(4, 0, 1.6);
  Vector3d camera_lookat(0, 0, 1);
  Vector3d camera_vertical(0, 0, 1);
  camera->set(cVector3d(camera_pos), cVector3d(camera_lookat),
              cVector3d(camera_vertical));

  graphics->SetBackgroundColor({0.3, 0.5, 0.7});

  // Parse model
  auto tree_parser = TreeParser(model_file);
  std::unique_ptr<TreeKinematic> tree(tree_parser.loadDescToTree());

  auto tree_visual = new spline_sim::TreeVisual(tree.get());
  graphics->AddOwning(tree_visual);

  auto branch_material = cMaterial::create();
  branch_material->m_diffuse = cColorf(0.77, 0.75, 0.62);
  branch_material->m_ambient = cColorf(0.04, 0.01, 0.01);
  branch_material->m_specular = cColorf(0.0, 0.05, 0.05);
  branch_material->setShininess(100);
  tree_visual->branchMaterialIs(branch_material);

  auto fruit_material = cMaterial::create();
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
  GLFWmonitor *primary = glfwGetPrimaryMonitor();
  const GLFWvidmode *mode = glfwGetVideoMode(primary);

  // information about computer screen and GLUT display window
  int screenW = mode->width;
  int screenH = mode->height;
  int windowW = 0.8 * screenH;
  int windowH = 0.5 * screenH;
  int windowPosY = (screenH - windowH) / 2;
  int windowPosX = windowPosY;

  // create window and make it current
  glfwWindowHint(GLFW_VISIBLE, 0);
  GLFWwindow *window =
      glfwCreateWindow(windowW, windowH, "04-cable", NULL, NULL);
  glfwSetWindowPos(window, windowPosX, windowPosY);
  glfwShowWindow(window);
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  // set callbacks
  glfwSetKeyCallback(window, keySelect);
  glfwSetMouseButtonCallback(window, mouseClick);

  std::thread update_thread(update, tree.get());

  /*------- Loop -------*/
  // cache variables
  double last_cursorx, last_cursory;

  Eigen::MatrixXd G;
  Eigen::Matrix3d R;
  Eigen::Vector3d center_point = Eigen::Vector3d::Zero();

  // while window is open:
  while (!glfwWindowShouldClose(window)) {
    // update graphics. this automatically waits for the correct amount of time
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    // render scene
    tree_visual->updateGraphics();
    graphics->UpdateShadowMaps(false);
    graphics->Render(kCameraName, width, height);

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

    Eigen::Vector3d cam_up_axis;
    // cam_up_axis = camera_vertical;
    // cam_up_axis.normalize();
    cam_up_axis << 0.0, 0.0, 1.0; // TODO: there might be a better way to do
                                  // this
    Eigen::Vector3d cam_roll_axis =
        (camera_lookat - camera_pos).cross(cam_up_axis);
    cam_roll_axis.normalize();
    Eigen::Vector3d cam_lookat_axis = camera_lookat;
    cam_lookat_axis.normalize();
    if (fTransXp) {
      camera_pos = camera_pos + 0.05 * cam_roll_axis;
      camera_lookat = camera_lookat + 0.05 * cam_roll_axis;
    }
    if (fTransXn) {
      camera_pos = camera_pos - 0.05 * cam_roll_axis;
      camera_lookat = camera_lookat - 0.05 * cam_roll_axis;
    }
    if (fTransYp) {
      // camera_pos = camera_pos + 0.05*cam_lookat_axis;
      camera_pos = camera_pos + 0.05 * cam_up_axis;
      camera_lookat = camera_lookat + 0.05 * cam_up_axis;
    }
    if (fTransYn) {
      // camera_pos = camera_pos - 0.05*cam_lookat_axis;
      camera_pos = camera_pos - 0.05 * cam_up_axis;
      camera_lookat = camera_lookat - 0.05 * cam_up_axis;
    }
    if (fRotPanTilt) {
      // get current cursor position
      double cursorx, cursory;
      glfwGetCursorPos(window, &cursorx, &cursory);
      // TODO: might need to re-scale from screen units to physical units
      double compass = 0.006 * (cursorx - last_cursorx);
      double azimuth = 0.006 * (cursory - last_cursory);
      double radius = (camera_pos - camera_lookat).norm();
      Eigen::Matrix3d m_tilt;
      m_tilt = Eigen::AngleAxisd(azimuth, -cam_roll_axis);
      camera_pos = camera_lookat + m_tilt * (camera_pos - camera_lookat);
      Eigen::Matrix3d m_pan;
      m_pan = Eigen::AngleAxisd(compass, -cam_up_axis);
      camera_pos = camera_lookat + m_pan * (camera_pos - camera_lookat);
    }
    camera->set(cVector3d(camera_pos), cVector3d(camera_lookat),
                cVector3d(camera_vertical));
    glfwGetCursorPos(window, &last_cursorx, &last_cursory);
  }

  // stop simulation
  fSimulationRunning = false;
  update_thread.join();

  // destroy context
  glfwDestroyWindow(window);

  // terminate
  glfwTerminate();

  return 0;
}

//------------------------------------------------------------------------------
void update(TreeKinematic *tree_kinematic) {
  // create a timer
  LoopTimer timer;
  timer.initializeTimer();
  timer.setLoopFrequency(1000);           // 1000Hz timer
  double last_time = timer.elapsedTime(); // secs

  bool fTimerDidSleep = true;

  spline_sim::TreeDynamics dynamics(tree_kinematic);

  // start simulation loop
  fSimulationRunning = true;
  while (fSimulationRunning) { // automatically set to false when simulation is
                               // quit
    fTimerDidSleep = timer.waitForNextLoop();

    // update time
    double curr_time = timer.elapsedTime();
    double loop_dt = curr_time - last_time;

    dynamics.Step(loop_dt, /*contact_list=*/{});

    // -------------------------------------------
    // update last time
    last_time = curr_time;
  }
}

//------------------------------------------------------------------------------

void glfwError(int error, const char *description) {
  std::cerr << "GLFW Error: " << description << std::endl;
  exit(1);
}

//------------------------------------------------------------------------------

void keySelect(GLFWwindow *window, int key, int scancode, int action,
               int mods) {
  bool set = (action != GLFW_RELEASE);
  switch (key) {
  case GLFW_KEY_ESCAPE:
    // exit application
    glfwSetWindowShouldClose(window, GL_TRUE);
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

void mouseClick(GLFWwindow *window, int button, int action, int mods) {
  bool set = (action != GLFW_RELEASE);
  // TODO: mouse interaction with robot
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
    // TODO: menu
    break;
  // if middle click: don't handle. doesn't work well on laptops
  case GLFW_MOUSE_BUTTON_MIDDLE:
    break;
  default:
    break;
  }
}
