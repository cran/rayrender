#include "disk.h"


bool disk::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 n(0.0, 1.0, 0.0);
  // First we intersect with the plane containing the disk
  Float t = -r.origin().y() / r.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float z = r.origin().z() + t*r.direction().z();
  Float radHit2 = x*x + z*z;
  if(radHit2 >= radius * radius || radHit2 <= inner_radius * inner_radius) {
    return(false);
  }
  
  
  vec3 p = r.point_at_parameter(t);
  p.e[1] = 0;
  
  Float u = p.x() / (2.0 * radius) + 0.5;
  Float v = p.z() / (2.0 * radius) + 0.5;
  u = 1 - u;
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < rng.unif_rand()) {
      return(false);
    }
  }
  rec.p = p;
  rec.normal = n;
  rec.t = t;
  rec.mat_ptr = mat_ptr.get();
  rec.u = u;
  rec.v = v;
  
  //Interaction information
  rec.dpdu = vec3(1, 0, 0);
  rec.dpdv = vec3(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  return(true);
}


bool disk::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  vec3 n(0.0, 1.0, 0.0);
  // First we intersect with the plane containing the disk
  Float t = -r.origin().y() / r.direction().y();
  if(t < t_min || t > t_max) {
    return(false);
  }
  Float x = r.origin().x() + t*r.direction().x();
  Float z = r.origin().z() + t*r.direction().z();
  Float radHit2 = x*x + z*z;
  if(radHit2 >= radius * radius || radHit2 <= inner_radius * inner_radius) {
    return(false);
  }
  
  
  vec3 p = r.point_at_parameter(t);
  p.e[1] = 0;
  
  Float u = p.x() / (2.0 * radius) + 0.5;
  Float v = p.z() / (2.0 * radius) + 0.5;
  u = 1 - u;
  if(alpha_mask) {
    if(alpha_mask->value(u, v, rec.p).x() < sampler->Get1D()) {
      return(false);
    }
  }
  rec.p = p;
  rec.normal = n;
  rec.t = t;
  rec.mat_ptr = mat_ptr.get();
  rec.u = u;
  rec.v = v;
  
  //Interaction information
  rec.dpdu = vec3(1, 0, 0);
  rec.dpdv = vec3(0, 0, 1);
  rec.has_bump = bump_tex ? true : false;
  
  if(bump_tex) {
    vec3 bvbu = bump_tex->value(rec.u,rec.v, rec.p);
    rec.bump_normal = rec.normal + bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv; 
    rec.bump_normal.make_unit_vector();
  }
  
  return(true);
}


Float disk::pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, rng)) {
    Float area =  M_PI * (radius * radius - inner_radius * inner_radius);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

Float disk::pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time) {
  hit_record rec;
  if(this->hit(ray(o,v), 0.001, FLT_MAX, rec, sampler)) {
    Float area =  M_PI * (radius * radius - inner_radius * inner_radius);
    Float distance_squared = rec.t * rec.t * v.squared_length();
    Float cosine = fabs(dot(v,rec.normal)/v.length());
    return(distance_squared / (cosine * area));
  } else {
    return(0);
  }
}

vec3 disk::random(const vec3& o, random_gen& rng, Float time) {
  Float r1 = rng.unif_rand();
  Float r2 = sqrt(rng.unif_rand());
  Float phi = 2 * M_PI * r1;
  Float x = ((radius - inner_radius) * r2 + inner_radius) * cos(phi);
  Float z = ((radius - inner_radius) * r2 + inner_radius) * sin(phi);
  return(vec3(x,0,z)+center-o);
}

vec3 disk::random(const vec3& o, Sampler* sampler, Float time) {
  vec2 u = sampler->Get2D();
  Float r1 = u.x();
  Float r2 = sqrt(u.y());
  Float phi = 2 * M_PI * r1;
  Float x = ((radius - inner_radius) * r2 + inner_radius) * cos(phi);
  Float z = ((radius - inner_radius) * r2 + inner_radius) * sin(phi);
  return(vec3(x,0,z)+center-o);
}

bool disk::bounding_box(Float t0, Float t1, aabb& box) const {
  box = aabb(-vec3(radius,0.001,radius), vec3(radius,0.001,radius));
  return(true);
}
