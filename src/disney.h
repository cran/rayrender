#ifndef TEXTUREH
#define TEXTUREH

#include "texture.h"
#include "vec3.h"
#include "mathinline.h"

inline Float AbsCosTheta(vec3 w) {
  return(std::fabs(w.z()));
}

class disney_texture : public texture {
public: 
  disney_texture() {}
  disney_texture(vec3 c, Float sheen_, vec3 sheenTint_) : baseColor(c), sheen(sheen_), sheenTint(sheenTint_) {
    Float luminance = dot(vec3(0.3f, 0.6f, 1.0f), baseColor);
    tint = (luminance > 0.0f) ? baseColor * (1.0f / luminance) : vec3(1,1,1);
  }
  virtual vec3 value(Float u, Float v, const vec3& p) const {
    return(baseColor);
  }
  const Float schlick(Float cosine) {
    Float r0 = (1 - 1.5) / (1 + 1.5);
    r0 = r0 * r0;
    return(r0 + (1-r0) * pow((1-cosine),5));
  }
  vec3 EvaluateSheen(const vec3& raydir, const vec3& wm, const vec3& wi) {
    if(sheen <= 0.0f) {
      return(vec3(0,0,0));
    }
    Float dotHL = dot(wm, wi);
    return(sheen * lerp(sheenTint, vec3(1.0f, 1.0f, 1.0f), tint) * schlick(dotHL));
  }
  static Float GTR1(Float absDotHL, Float a) {
    if(a >= 1) {
      return M_1_PI;
    }
    Float a2 = a * a;
    return((a2 - 1.0f) / (M_PI * log2(a2) * (1.0f + (a2 - 1.0f) * absDotHL * absDotHL)));
  }
  float SeparableSmithGGXG1(const vec3& w, Float a) {
    Float a2 = a * a;
    Float absDotNV = AbsCosTheta(w);
    
    return((2.0f / (1.0f + std::sqrt(a2 + (1 - a2) * absDotNV * absDotNV))));
  }
  
  vec3 baseColor, tint, sheenTint;
  Float sheen;
};

#endif
