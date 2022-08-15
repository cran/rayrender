#ifndef SINGLETRIANGLEH
#define SINGLETRIANGLEH

#include "trianglemesh.h"
#include "transform.h"
#include "hitable.h"

class Triangle : public hitable {
public:
Triangle(std::shared_ptr<Transform>ObjectToWorld, std::shared_ptr<Transform> WorldToObject,
         bool reverseOrientation,
         const std::shared_ptr<TriangleMesh> &mesh, int triNumber)
  : hitable(ObjectToWorld, WorldToObject, reverseOrientation),
    mesh(mesh) {
  v = &mesh->vertexIndices[3 * triNumber];
}
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler);
  
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
  virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time = 0);
  virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time = 0);
  
  virtual vec3f random(const point3f& origin, random_gen& rng, Float time = 0);
  virtual vec3f random(const point3f& origin, Sampler* sampler, Float time = 0);
  virtual std::string GetName() const {
    return(std::string("Triangle"));
  }
  Float Area() const;
  
  std::shared_ptr<TriangleMesh> mesh;
  const int *v;
  std::shared_ptr<material> mp;
  std::shared_ptr<bump_texture> bump_tex;
  std::shared_ptr<alpha_texture> alpha_mask;
  
private:
  void GetUVs(point2f uv[3]) const {
    if (mesh->uv) {
      uv[0] = mesh->uv[v[0]];
      uv[1] = mesh->uv[v[1]];
      uv[2] = mesh->uv[v[2]];
    } else {
      uv[0] = point2f(0, 0);
      uv[1] = point2f(1, 0);
      uv[2] = point2f(1, 1);
    }
}
  
  
};

#endif 
