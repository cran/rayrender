#include "hitable.h"

//Translate implementation

void get_sphere_uv(const vec3f& p, Float& u, Float& v) {
  Float phi = atan2(p.z(),p.x());
  Float theta = asin(p.y());
  u = 1 - (phi + M_PI) / (2*M_PI);
  v = (theta + M_PI/2) / M_PI;
}

void get_sphere_uv(const normal3f& p, Float& u, Float& v) {
  Float phi = atan2(p.z(),p.x());
  Float theta = asin(p.y());
  u = 1 - (phi + M_PI) / (2*M_PI);
  v = (theta + M_PI/2) / M_PI;
}

Float AnimatedHitable::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  return(primitive->pdf_value(Inverse(InterpolatedPrimToWorld)(o),
                              Inverse(InterpolatedPrimToWorld)(v),
                              rng, time));
}
Float AnimatedHitable::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  return(primitive->pdf_value(Inverse(InterpolatedPrimToWorld)(o),
                              Inverse(InterpolatedPrimToWorld)(v),
                              sampler, time));
}
vec3f AnimatedHitable::random(const point3f& o, random_gen& rng, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  return(InterpolatedPrimToWorld(primitive->random(Inverse(InterpolatedPrimToWorld)(o), rng, time)));
}
vec3f AnimatedHitable::random(const point3f& o, Sampler* sampler, Float time) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(time, &InterpolatedPrimToWorld);
  
  return(InterpolatedPrimToWorld(primitive->random(Inverse(InterpolatedPrimToWorld)(o), sampler, time)));
}
std::string AnimatedHitable::GetName() const {
  return(std::string("AnimatedHitable"));
}



bool AnimatedHitable::bounding_box(Float t0, Float t1, aabb& box) const {
  primitive->bounding_box(t0, t1, box);
  box = PrimitiveToWorld.MotionBounds(box);
  return(true);
}

bool AnimatedHitable::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(r.time(), &InterpolatedPrimToWorld);
  ray ray_interp = Inverse(InterpolatedPrimToWorld)(r);
  
  if (!primitive->hit(ray_interp, t_min, t_max, rec, rng)) {
    return false;
  }
  r.tMax = ray_interp.tMax;
  if (!InterpolatedPrimToWorld.IsIdentity()) {
    rec = InterpolatedPrimToWorld(rec);
  }
  return true;
}

bool AnimatedHitable::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  Transform InterpolatedPrimToWorld;
  PrimitiveToWorld.Interpolate(r.time(), &InterpolatedPrimToWorld);
  ray ray_interp = Inverse(InterpolatedPrimToWorld)(r);
  
  if (!primitive->hit(ray_interp, t_min, t_max, rec, sampler)) {
    return false;
  }
  r.tMax = ray_interp.tMax;
  if (!InterpolatedPrimToWorld.IsIdentity()) {
    rec = InterpolatedPrimToWorld(rec);
  }
  return true;
}
