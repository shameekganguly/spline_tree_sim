#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <memory>

#include <chai3d.h>

namespace spline_sim {

class Graphics {
public:
  Graphics() : chai3d_world_(std::make_unique<chai3d::cWorld>()) {}

  bool AddOwning(chai3d::cGenericObject *object);

  chai3d::cCamera *CreateCamera(std::string_view name);

  void CreateLight(const chai3d::cVector3d &pos,
                   const chai3d::cVector3d &look_at);

  void SetBackgroundColor(const std::vector<double> &rgb);

  void UpdateShadowMaps(bool update);

  void Render(std::string_view camera_name, int width, int height);

protected:
  std::unique_ptr<chai3d::cWorld> chai3d_world_;
  std::unordered_map<std::string, chai3d::cCamera *> cameras_;
};
} // namespace spline_sim

#endif // GRAPHICS_H
