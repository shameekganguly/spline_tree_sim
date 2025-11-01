#include "Graphics.h"

namespace spline_sim {

bool Graphics::AddOwning(chai3d::cGenericObject *object) {
  chai3d_world_->addChild(object);
  return true;
}

chai3d::cCamera *Graphics::CreateCamera(std::string_view name) {
  assert(chai3d_world_ != nullptr);
  chai3d::cCamera *camera = new chai3d::cCamera(chai3d_world_.get());
  camera->m_name = name;
  // Transfer ownership
  chai3d_world_->addChild(camera);

  // Default properties
  camera->setClippingPlanes(0.01, 10.0);

  cameras_[std::string(name)] = camera;
  return camera;
}

void Graphics::SetBackgroundColor(const std::vector<double> &rgb) {
  assert(rgb.size() == 3);
  chai3d_world_->setBackgroundColor(rgb[0], rgb[1], rgb[2]);
}

void Graphics::UpdateShadowMaps(bool update) {
  chai3d_world_->updateShadowMaps(update);
}

void Graphics::Render(std::string_view camera_name, int width, int height) {
  cameras_[std::string(camera_name)]->renderView(width, height);
}

} // namespace spline_sim
