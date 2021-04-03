#include "csg.h"

bool csg::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  Float threshold = 0.001;
  
  Float delta = 10e-5 * max_dist/100; 
  Float t = 0; 
  bool first = true;
  aabb box;
  bounding_box(t_min,t_max,box);
  Float dist1 = (box.min()-r.origin()).length();
  Float dist2 = (box.max()-r.origin()).length();
  
  Float max_path = std::fmax(dist1,dist2) + box.diag.length()/2;
  
  vec3 dir = unit_vector(r.direction());
  Float max_t = max_path * r.direction().length();
  
  
  while (t < max_t) { 
    Float minDistance = INFINITY; 
    vec3 from = r.origin() + t * dir; 
    
    //Need distance from interior edge to support dielectrics
    float d =  std::fabs(shapes->getDistance(from)); 
    
    //Need to deal with refraction, often initial distance is too close to surface, so we offset
    if(first && d < threshold) {
      t += 0.01; //Hard coded offset, not great
      first = false;
      continue;
    } 
    
    if (d < minDistance) {
      minDistance = d;
    }
    first = false;
    if (minDistance <= threshold) { 
      Float tval = t / r.direction().length();
      if(tval > t_min && tval < t_max) {
        rec.normal = vec3( 
          shapes->getDistance(from + vec3(delta, 0, 0)) - shapes->getDistance(from + vec3(-delta, 0, 0)), 
          shapes->getDistance(from + vec3(0, delta, 0)) - shapes->getDistance(from + vec3(0, -delta, 0)), 
          shapes->getDistance(from + vec3(0, 0, delta)) - shapes->getDistance(from + vec3(0, 0, -delta))
        );
        //Deal with degenerate case by setting directly at camera--not ideal, need better fix
        if(rec.normal.x() == 0 && rec.normal.y() == 0 && rec.normal.z() == 0) {
          rec.normal = -r.direction();
        }
        rec.p = from;
        rec.normal.make_unit_vector(); 
        rec.t = tval;
        rec.u = 0.5;
        rec.v = 0.5;
        rec.dpdu = 0.5;
        rec.dpdv = 0.5;
        rec.mat_ptr = mat_ptr.get();
        rec.has_bump = false;
        rec.bump_normal =  rec.normal;
        return(true);
      } else {
        return(false);
      }
    } 
    t += minDistance; 
  } 
  return(false);
}


bool csg::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  Float threshold = 0.001;
  
  Float delta = 10e-5 * max_dist/100; 
  Float t = 0; 
  bool first = true;
  vec3 dir = unit_vector(r.direction());
  Float max_t = t_max * r.direction().length();
  
  
  while (t < max_t) { 
    Float minDistance = INFINITY; 
    vec3 from = r.origin() + t * dir; 
    
    //Need distance from interior edge to support dielectrics
    float d =  std::fabs(shapes->getDistance(from)); 
    
    //Need to deal with refraction, often initial distance is too close to surface, so we offset
    if(first && d < threshold) {
      t += 0.01; //Hard coded offset, not great
      first = false;
      continue;
    } 
    
    if (d < minDistance) {
      minDistance = d;
    }
    first = false;
    if (minDistance <= threshold) { 
      Float tval = t / r.direction().length();
      if(tval > t_min && tval < t_max) {
        rec.normal = vec3( 
          shapes->getDistance(from + vec3(delta, 0, 0)) - shapes->getDistance(from + vec3(-delta, 0, 0)), 
          shapes->getDistance(from + vec3(0, delta, 0)) - shapes->getDistance(from + vec3(0, -delta, 0)), 
          shapes->getDistance(from + vec3(0, 0, delta)) - shapes->getDistance(from + vec3(0, 0, -delta))
        );
        //Deal with degenerate case by setting directly at camera--not ideal, need better fix
        if(rec.normal.x() == 0 && rec.normal.y() == 0 && rec.normal.z() == 0) {
          rec.normal = -r.direction();
        }
        rec.p = from;
        rec.normal.make_unit_vector(); 
        rec.t = tval;
        rec.u = 0.5;
        rec.v = 0.5;
        rec.dpdu = 0.5;
        rec.dpdv = 0.5;
        rec.mat_ptr = mat_ptr.get();
        rec.has_bump = false;
        rec.bump_normal =  rec.normal;
        return(true);
      } else {
        return(false);
      }
    } 
    t += minDistance; 
  } 
  return(false);
}

bool csg::bounding_box(Float t0, Float t1, aabb& box) const {
  shapes->bbox(t0,t1,box);
  return(true);
}
