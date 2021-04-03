#include "triangle.h"


bool triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 pvec = cross(r.direction(), edge2);
  Float det = dot(pvec, edge1);
  
  // no culling
  if (std::fabs(det) < 1E-15) {
    return(false);
  }
  Float invdet = 1.0 / det;
  vec3 tvec = r.origin() - a;
  Float u = dot(pvec, tvec) * invdet;
  if (u < 0.0 || u > 1.0) {
    return(false);
  }
  
  vec3 qvec = cross(tvec, edge1);
  Float v = dot(qvec, r.direction()) * invdet;
  if (v < 0 || u + v > 1.0) {
    return(false);
  }
  Float t = dot(qvec, edge2) * invdet; 
  
  if (t < t_min || t > t_max) {
    return(false);
  }
  if(alpha_mask) {
    if(alpha_mask->channel_value(u, v, rec.p) < rng.unif_rand()) {
      return(false);
    }
  }
  Float w = 1 - u - v;
  rec.t = t;
  rec.p = r.point_at_parameter(t);
  rec.u = u;
  rec.v = v;
  rec.has_bump = false;
  
  if(bump_tex) {
    //Get UV values + calculate dpdu/dpdv
    vec3 u_val, v_val;
    vec2 uv[3];
    u_val = bump_tex->u_vec;
    v_val = bump_tex->v_vec;
    uv[0] = vec2(u_val.x(), v_val.x());
    uv[1] = vec2(u_val.y(), v_val.y());
    uv[2] = vec2(u_val.z(), v_val.z());
    vec2 duv02 = uv[0] - uv[2];
    vec2 duv12 = uv[1] - uv[2];
    vec3 dp02 = edge1;
    vec3 dp12 = edge2;
    
    Float determinant = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
    if (determinant == 0) {
      onb uvw;
      uvw.build_from_w(cross(edge2,edge1));
      rec.dpdu = uvw.u();
      rec.dpdv = uvw.v();
    } else {
      Float invdet = 1 / determinant;
      rec.dpdu = -( duv12[1] * dp02 - duv02[1] * dp12) * invdet;
      rec.dpdv = -(-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
  } else {
    rec.dpdu = vec3(0,0,0);
    rec.dpdv = vec3(0,0,0);
  }
  // Use that to calculate normals
  if(normals_provided) {
    vec3 normal_temp = w * na + u * nb + v * nc;
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal_temp) < 0 ? normal_temp : -normal_temp;
    } else {
      rec.normal = normal_temp;
    }
  } else {
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal) < 0 ? normal : -normal;
    } else {
      rec.normal = normal;
    }
  }
  if(bump_tex) {
    vec3 bvbu = bump_tex->mesh_value(rec.u, rec.v, rec.p);
    rec.bump_normal = dot(r.direction(), normal) < 0 ? 
    rec.normal +  bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv :
      rec.normal -  bvbu.x() * rec.dpdu - bvbu.y() * rec.dpdv;
    rec.bump_normal.make_unit_vector();
    rec.has_bump = true;
  }
  rec.mat_ptr = mp.get();
  return(true);
}


bool triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  vec3 pvec = cross(r.direction(), edge2);
  Float det = dot(pvec, edge1);
  
  // no culling
  if (std::fabs(det) < 1E-15) {
    return(false);
  }
  Float invdet = 1.0 / det;
  vec3 tvec = r.origin() - a;
  Float u = dot(pvec, tvec) * invdet;
  if (u < 0.0 || u > 1.0) {
    return(false);
  }
  
  vec3 qvec = cross(tvec, edge1);
  Float v = dot(qvec, r.direction()) * invdet;
  if (v < 0 || u + v > 1.0) {
    return(false);
  }
  Float t = dot(qvec, edge2) * invdet; 
  
  if (t < t_min || t > t_max) {
    return(false);
  }
  if(alpha_mask) {
    if(alpha_mask->channel_value(u, v, rec.p) < sampler->Get1D()) {
      return(false);
    }
  }
  Float w = 1 - u - v;
  rec.t = t;
  rec.p = r.point_at_parameter(t);
  rec.u = u;
  rec.v = v;
  rec.has_bump = false;
  
  if(bump_tex) {
    //Get UV values + calculate dpdu/dpdv
    vec3 u_val, v_val;
    vec2 uv[3];
    u_val = bump_tex->u_vec;
    v_val = bump_tex->v_vec;
    uv[0] = vec2(u_val.x(), v_val.x());
    uv[1] = vec2(u_val.y(), v_val.y());
    uv[2] = vec2(u_val.z(), v_val.z());
    vec2 duv02 = uv[0] - uv[2];
    vec2 duv12 = uv[1] - uv[2];
    vec3 dp02 = edge1;
    vec3 dp12 = edge2;
    
    Float determinant = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
    if (determinant == 0) {
      onb uvw;
      uvw.build_from_w(cross(edge2,edge1));
      rec.dpdu = uvw.u();
      rec.dpdv = uvw.v();
    } else {
      Float invdet = 1 / determinant;
      rec.dpdu = -( duv12[1] * dp02 - duv02[1] * dp12) * invdet;
      rec.dpdv = -(-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
  } else {
    rec.dpdu = vec3(0,0,0);
    rec.dpdv = vec3(0,0,0);
  }
  // Use that to calculate normals
  if(normals_provided) {
    vec3 normal_temp = w * na + u * nb + v * nc;
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal_temp) < 0 ? normal_temp : -normal_temp;
    } else {
      rec.normal = normal_temp;
    }
  } else {
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal) < 0 ? normal : -normal;
    } else {
      rec.normal = normal;
    }
  }
  if(bump_tex) {
    vec3 bvbu = bump_tex->mesh_value(rec.u, rec.v, rec.p);
    rec.bump_normal = dot(r.direction(), normal) < 0 ? 
    rec.normal +  bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv :
      rec.normal -  bvbu.x() * rec.dpdu - bvbu.y() * rec.dpdv;
    rec.bump_normal.make_unit_vector();
    rec.has_bump = true;
  }
  rec.mat_ptr = mp.get();
  return(true);
}

bool triangle::bounding_box(Float t0, Float t1, aabb& box) const {
  vec3 min_v(fmin(fmin(a.x(), b.x()), c.x()), 
             fmin(fmin(a.y(), b.y()), c.y()), 
             fmin(fmin(a.z(), b.z()), c.z()));
  vec3 max_v(fmax(fmax(a.x(), b.x()), c.x()), 
             fmax(fmax(a.y(), b.y()), c.y()), 
             fmax(fmax(a.z(), b.z()), c.z()));
  
  vec3 difference = max_v - min_v;
  
  if (difference.x() < 1E-5) max_v.e[0] += 1E-5;
  if (difference.y() < 1E-5) max_v.e[1] += 1E-5;
  if (difference.z() < 1E-5) max_v.e[2] += 1E-5;
  
  box = aabb(min_v, max_v);
  return(true);
}

Float triangle::pdf_value(const vec3& o, const vec3& v, random_gen& rng, Float time) { 
  hit_record rec;
  if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, rng)) {
    Float distance = rec.t * rec.t * v.squared_length();;
    Float cosine = dot(v, rec.normal);
    return(distance / (cosine * area));
  }
  return 0; 
}

Float triangle::pdf_value(const vec3& o, const vec3& v, Sampler* sampler, Float time) { 
  hit_record rec;
  if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, sampler)) {
    Float distance = rec.t * rec.t * v.squared_length();;
    Float cosine = dot(v, rec.normal);
    return(distance / (cosine * area));
  }
  return 0; 
}

vec3 triangle::random(const vec3& origin, random_gen& rng, Float time) {
  Float r1 = rng.unif_rand();
  Float r2 = rng.unif_rand();
  Float sr1 = sqrt(r1);
  vec3 random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin); 
}
vec3 triangle::random(const vec3& origin, Sampler* sampler, Float time) {
  vec2 u = sampler->Get2D();
  Float r1 = u.x();
  Float r2 = u.y();
  Float sr1 = sqrt(r1);
  vec3 random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin); 
}
