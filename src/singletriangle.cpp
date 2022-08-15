#include "singletriangle.h"


bool Triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  // Get Triangle vertices in _p0_, _p1_, and _p2_
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  // Perform ray--Triangle intersection test
  
  // Transform Triangle vertices to ray coordinate space
  
  // Translate vertices based on ray origin
  point3f p0t = p0 - vec3f(r.origin());
  point3f p1t = p1 - vec3f(r.origin());
  point3f p2t = p2 - vec3f(r.origin());
  
  // Permute components of Triangle vertices and ray direction
  int kz = MaxDimension(Abs(r.direction()));
  int kx = kz + 1;
  if (kx == 3) kx = 0;
  int ky = kx + 1;
  if (ky == 3) ky = 0;
  vec3f d = Permute(r.direction(), kx, ky, kz);
  p0t = Permute(p0t, kx, ky, kz);
  p1t = Permute(p1t, kx, ky, kz);
  p2t = Permute(p2t, kx, ky, kz);
  
  // Apply shear transformation to translated vertex positions
  Float Sx = -d.x() / d.z();
  Float Sy = -d.y() / d.z();
  Float Sz = 1.f /  d.z();
  p0t.e[0] += Sx * p0t.z();
  p0t.e[1] += Sy * p0t.z();
  p1t.e[0] += Sx * p1t.z();
  p1t.e[1] += Sy * p1t.z();
  p2t.e[0] += Sx * p2t.z();
  p2t.e[1] += Sy * p2t.z();
  
  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = p1t.x() * p2t.y() - p1t.y() * p2t.x();
  Float e1 = p2t.x() * p0t.y() - p2t.y() * p0t.x();
  Float e2 = p0t.x() * p1t.y() - p0t.y() * p1t.x();
  
  // Fall back to double precision test at Triangle edges
  if (sizeof(Float) == sizeof(float) &&
      (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
    double p2txp1ty = (double)p2t.x() * (double)p1t.y();
    double p2typ1tx = (double)p2t.y() * (double)p1t.x();
    e0 = (float)(p2typ1tx - p2txp1ty);
    double p0txp2ty = (double)p0t.x() * (double)p2t.y();
    double p0typ2tx = (double)p0t.y() * (double)p2t.x();
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.x() * (double)p0t.y();
    double p1typ0tx = (double)p1t.y() * (double)p0t.x();
    e2 = (float)(p1typ0tx - p1txp0ty);
  }
  
  // Perform Triangle edge and determinant tests
  if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
    return false;
  Float det = e0 + e1 + e2;
  if (det == 0) return false;
  
  // Compute scaled hit distance to Triangle and test against ray $t$ range
  p0t.e[2] *= Sz;
  p1t.e[2] *= Sz;
  p2t.e[2] *= Sz;
  Float tScaled = e0 * p0t.z() + e1 * p1t.z() + e2 * p2t.z();
  if (det < 0 && (tScaled >= 0 || tScaled < t_max * det)) {
    return false;
  }
  else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det)) {
    return false;
  }
  
  // Compute barycentric coordinates and $t$ value for Triangle intersection
  Float invDet = 1 / det;
  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  Float t = tScaled * invDet;
  
  // Ensure that computed Triangle $t$ is conservatively greater than zero
  
  // Compute $\delta_z$ term for Triangle $t$ error bounds
  Float maxZt = MaxComponent(Abs(vec3f(p0t.z(), p1t.z(), p2t.z())));
  Float deltaZ = gamma(3) * maxZt;
  
  // Compute $\delta_x$ and $\delta_y$ terms for Triangle $t$ error bounds
  Float maxXt = MaxComponent(Abs(vec3f(p0t.x(), p1t.x(), p2t.x())));
  Float maxYt = MaxComponent(Abs(vec3f(p0t.y(), p1t.y(), p2t.y())));
  Float deltaX = gamma(5) * (maxXt + maxZt);
  Float deltaY = gamma(5) * (maxYt + maxZt);
  
  // Compute $\delta_e$ term for Triangle $t$ error bounds
  Float deltaE =
    2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);
  
  // Compute $\delta_t$ term for Triangle $t$ error bounds and check _t_
  Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
  Float deltaT = 3 *
    (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
    std::abs(invDet);
  if (t <= deltaT) return false;
  
  // Compute Triangle partial derivatives
  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);
  
  // Compute deltas for Triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
  bool degenerateUV = std::abs(determinant) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinant;
    dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || cross(dpdu, dpdv).squared_length() == 0) {
    // Handle zero determinant for Triangle partial derivative matrix
    vec3f ng = cross(p2 - p0, p1 - p0);
    if (ng.squared_length() == 0)
      // The Triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    
    CoordinateSystem(unit_vector(ng), &dpdu, &dpdv);
  }
  rec.dpdu = dpdu;
  rec.dpdv = dpdv;
  
  // Compute error bounds for Triangle intersection
  Float xAbsSum =
    (std::abs(b0 * p0.x()) + std::abs(b1 * p1.x()) + std::abs(b2 * p2.x()));
  Float yAbsSum =
    (std::abs(b0 * p0.y()) + std::abs(b1 * p1.y()) + std::abs(b2 * p2.y()));
  Float zAbsSum =
    (std::abs(b0 * p0.z()) + std::abs(b1 * p1.z()) + std::abs(b2 * p2.z()));
  // vec3f pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  
  // Interpolate $(u,v)$ parametric coordinates and hit point
  point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
  
  Float u1 = uvHit[0];
  Float v1 = uvHit[1];
  // Test intersection against alpha texture, if present
  // if (testAlphaTexture && mesh->alphaMask) {
  //   SurfaceInteraction isectLocal(pHit, vec3f(0, 0, 0), uvHit, -r.direction(),
  //                                 dpdu, dpdv, normal3f(0, 0, 0),
  //                                 normal3f(0, 0, 0), ray.time, this);
  //   if (mesh->alphaMask->Evaluate(isectLocal) == 0) return false;
  // }
  bool alpha_miss = false;
  if(alpha_mask) {
    if(alpha_mask->value(u1, v1, rec.p) < rng.unif_rand()) {
      alpha_miss = true;
    }
  }
  
  // Fill in _SurfaceInteraction_ from Triangle hit
  // *isect = SurfaceInteraction(pHit, pError, uvHit, -r.direction(), dpdu, dpdv,
  //                             normal3f(0, 0, 0), normal3f(0, 0, 0), ray.time,
  //                             this, faceIndex);
  rec.t = t;
  rec.p = pHit;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  rec.has_bump = false;
  rec.u = 1-u1;
  rec.v = v1;
  rec.mat_ptr = mp.get();
  rec.alpha_miss = alpha_miss;
  rec.shape = this;
  
  // Override surface normal in _isect_ for Triangle
  // isect->n = isect->shading.n = normal3f(unit_vector(cross(dp02, dp12)));
  rec.normal = normal3f(unit_vector(cross(dp02, dp12)));
  
  if (reverseOrientation ^ transformSwapsHandedness) {
    rec.normal = -rec.normal;
  }
  
  if (mesh->n || mesh->s) {
    // Initialize _Triangle_ shading geometry
    
    // Compute shading normal _ns_ for Triangle
    normal3f ns;
    if (mesh->n) {
      ns = (b0 * mesh->n[v[0]] + b1 * mesh->n[v[1]] + b2 * mesh->n[v[2]]);
      if (ns.squared_length() > 0) {
          ns = unit_vector(ns);
      } else {
        ns = rec.normal;
      }
    } else {
      ns = rec.normal;
    }
    
    // Compute shading tangent _ss_ for Triangle
    vec3f ss;
    if (mesh->s) {
      ss = (b0 * mesh->s[v[0]] + b1 * mesh->s[v[1]] + b2 * mesh->s[v[2]]);
      if (ss.squared_length() > 0) {
        ss = unit_vector(ss);
      }  else {
        ss = unit_vector(rec.dpdu);
      }
    } else {
      ss = unit_vector(rec.dpdu);
    }
    
    // Compute shading bitangent _ts_ for Triangle and adjust _ss_
    vec3f ts = cross(ss, ns.convert_to_vec3());
    if (ts.squared_length() > 0.f) {
      ts = unit_vector(ts);
      ss = cross(ts, ns.convert_to_vec3());
    } else {
      CoordinateSystem(ns.convert_to_vec3(), &ss, &ts);
    }
    
    // Compute $\dndu$ and $\dndv$ for Triangle shading geometry
    normal3f dndu, dndv;
    if (mesh->n) {
      // Compute deltas for Triangle partial derivatives of normal
      vec2f duv02 = uv[0] - uv[2];
      vec2f duv12 = uv[1] - uv[2];
      normal3f dn1 = mesh->n[v[0]] - mesh->n[v[2]];
      normal3f dn2 = mesh->n[v[1]] - mesh->n[v[2]];
      Float determinant = DifferenceOfProducts(duv02[0], duv12[1], duv02[1], duv12[0]);
      bool degenerateUV = std::abs(determinant) < 1e-8;
      if (degenerateUV) {
        // We can still compute dndu and dndv, with respect to the
        // same arbitrary coordinate system we use to compute dpdu
        // and dpdv when this happens. It's important to do this
        // (rather than giving up) so that ray differentials for
        // rays reflected from Triangles with degenerate
        // parameterizations are still reasonable.
        vec3f dn = cross((mesh->n[v[2]] - mesh->n[v[0]]).convert_to_vec3(),
                         (mesh->n[v[1]] - mesh->n[v[0]]).convert_to_vec3());
        if (dn.squared_length() == 0) {
          dndu = dndv = normal3f(0, 0, 0);
        } else {
          vec3f dnu, dnv;
          CoordinateSystem(dn, &dnu, &dnv);
          dndu = normal3f(dnu);
          dndv = normal3f(dnv);
        }
      } else {
        Float invDet = 1 / determinant;
        dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
        dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
      }
    } else {
      dndu = dndv = normal3f(0, 0, 0);
    }
    if (reverseOrientation) {
      ts = -ts;
    } 
    rec.dpdu = ss;
    rec.dpdv = ts;
    rec.dndu = dndu;
    rec.dndv = dndv;
    
    // isect->SetShadingGeometry(ss, ts, dndu, dndv, true);
  }
  rec.normal = unit_vector(cross(rec.dpdu, rec.dpdv));
  
  return true;
}

bool Triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  // Get Triangle vertices in _p0_, _p1_, and _p2_
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  // Perform ray--Triangle intersection test
  
  // Transform Triangle vertices to ray coordinate space
  
  // Translate vertices based on ray origin
  point3f p0t = p0 - vec3f(r.origin());
  point3f p1t = p1 - vec3f(r.origin());
  point3f p2t = p2 - vec3f(r.origin());
  
  // Permute components of Triangle vertices and ray direction
  int kz = MaxDimension(Abs(r.direction()));
  int kx = kz + 1;
  if (kx == 3) kx = 0;
  int ky = kx + 1;
  if (ky == 3) ky = 0;
  vec3f d = Permute(r.direction(), kx, ky, kz);
  p0t = Permute(p0t, kx, ky, kz);
  p1t = Permute(p1t, kx, ky, kz);
  p2t = Permute(p2t, kx, ky, kz);
  
  // Apply shear transformation to translated vertex positions
  Float Sx = -d.x() / d.z();
  Float Sy = -d.y() / d.z();
  Float Sz = 1.f /  d.z();
  p0t.e[0] += Sx * p0t.z();
  p0t.e[1] += Sy * p0t.z();
  p1t.e[0] += Sx * p1t.z();
  p1t.e[1] += Sy * p1t.z();
  p2t.e[0] += Sx * p2t.z();
  p2t.e[1] += Sy * p2t.z();
  
  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = p1t.x() * p2t.y() - p1t.y() * p2t.x();
  Float e1 = p2t.x() * p0t.y() - p2t.y() * p0t.x();
  Float e2 = p0t.x() * p1t.y() - p0t.y() * p1t.x();
  
  // Fall back to double precision test at Triangle edges
  if (sizeof(Float) == sizeof(float) &&
      (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
    double p2txp1ty = (double)p2t.x() * (double)p1t.y();
    double p2typ1tx = (double)p2t.y() * (double)p1t.x();
    e0 = (float)(p2typ1tx - p2txp1ty);
    double p0txp2ty = (double)p0t.x() * (double)p2t.y();
    double p0typ2tx = (double)p0t.y() * (double)p2t.x();
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.x() * (double)p0t.y();
    double p1typ0tx = (double)p1t.y() * (double)p0t.x();
    e2 = (float)(p1typ0tx - p1txp0ty);
  }
  
  // Perform Triangle edge and determinant tests
  if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
    return false;
  Float det = e0 + e1 + e2;
  if (det == 0) return false;
  
  // Compute scaled hit distance to Triangle and test against ray $t$ range
  p0t.e[2] *= Sz;
  p1t.e[2] *= Sz;
  p2t.e[2] *= Sz;
  Float tScaled = e0 * p0t.z() + e1 * p1t.z() + e2 * p2t.z();
  if (det < 0 && (tScaled >= 0 || tScaled < t_max * det)) {
    return false;
  }
  else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det)) {
    return false;
  }
  
  // Compute barycentric coordinates and $t$ value for Triangle intersection
  Float invDet = 1 / det;
  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  Float t = tScaled * invDet;
  
  // Ensure that computed Triangle $t$ is conservatively greater than zero
  
  // Compute $\delta_z$ term for Triangle $t$ error bounds
  Float maxZt = MaxComponent(Abs(vec3f(p0t.z(), p1t.z(), p2t.z())));
  Float deltaZ = gamma(3) * maxZt;
  
  // Compute $\delta_x$ and $\delta_y$ terms for Triangle $t$ error bounds
  Float maxXt = MaxComponent(Abs(vec3f(p0t.x(), p1t.x(), p2t.x())));
  Float maxYt = MaxComponent(Abs(vec3f(p0t.y(), p1t.y(), p2t.y())));
  Float deltaX = gamma(5) * (maxXt + maxZt);
  Float deltaY = gamma(5) * (maxYt + maxZt);
  
  // Compute $\delta_e$ term for Triangle $t$ error bounds
  Float deltaE =
    2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);
  
  // Compute $\delta_t$ term for Triangle $t$ error bounds and check _t_
  Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
  Float deltaT = 3 *
    (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
    std::abs(invDet);
  if (t <= deltaT) return false;
  
  // Compute Triangle partial derivatives
  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);
  
  // Compute deltas for Triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
  bool degenerateUV = std::abs(determinant) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinant;
    dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || cross(dpdu, dpdv).squared_length() == 0) {
    // Handle zero determinant for Triangle partial derivative matrix
    vec3f ng = cross(p2 - p0, p1 - p0);
    if (ng.squared_length() == 0)
      // The Triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    
    CoordinateSystem(unit_vector(ng), &dpdu, &dpdv);
  }
  rec.dpdu = dpdu;
  rec.dpdv = dpdv;
  
  // Compute error bounds for Triangle intersection
  Float xAbsSum =
    (std::abs(b0 * p0.x()) + std::abs(b1 * p1.x()) + std::abs(b2 * p2.x()));
  Float yAbsSum =
    (std::abs(b0 * p0.y()) + std::abs(b1 * p1.y()) + std::abs(b2 * p2.y()));
  Float zAbsSum =
    (std::abs(b0 * p0.z()) + std::abs(b1 * p1.z()) + std::abs(b2 * p2.z()));
  // vec3f pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  
  // Interpolate $(u,v)$ parametric coordinates and hit point
  point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
  
  Float u1 = uvHit[0];
  Float v1 = uvHit[1];
  // Test intersection against alpha texture, if present
  // if (testAlphaTexture && mesh->alphaMask) {
  //   SurfaceInteraction isectLocal(pHit, vec3f(0, 0, 0), uvHit, -r.direction(),
  //                                 dpdu, dpdv, normal3f(0, 0, 0),
  //                                 normal3f(0, 0, 0), ray.time, this);
  //   if (mesh->alphaMask->Evaluate(isectLocal) == 0) return false;
  // }
  bool alpha_miss = false;
  if(alpha_mask) {
    if(alpha_mask->value(u1, v1, rec.p) < sampler->Get1D()) {
      alpha_miss = true;
    }
  }
  
  // Fill in _SurfaceInteraction_ from Triangle hit
  // *isect = SurfaceInteraction(pHit, pError, uvHit, -r.direction(), dpdu, dpdv,
  //                             normal3f(0, 0, 0), normal3f(0, 0, 0), ray.time,
  //                             this, faceIndex);
  rec.t = t;
  rec.p = pHit;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  rec.has_bump = false;
  rec.u = 1-u1;
  rec.v = v1;
  rec.mat_ptr = mp.get();
  rec.alpha_miss = alpha_miss;
  rec.shape = this;
  
  // Override surface normal in _isect_ for Triangle
  // isect->n = isect->shading.n = normal3f(unit_vector(cross(dp02, dp12)));
  rec.normal = normal3f(unit_vector(cross(dp02, dp12)));
  
  if (reverseOrientation ^ transformSwapsHandedness) {
    rec.normal = -rec.normal;
  }
  
  if (mesh->n || mesh->s) {
    // Initialize _Triangle_ shading geometry
    
    // Compute shading normal _ns_ for Triangle
    normal3f ns;
    if (mesh->n) {
      ns = (b0 * mesh->n[v[0]] + b1 * mesh->n[v[1]] + b2 * mesh->n[v[2]]);
      if (ns.squared_length() > 0) {
        ns = unit_vector(ns);
      } else {
        ns = rec.normal;
      }
    } else {
      ns = rec.normal;
    }
    
    // Compute shading tangent _ss_ for Triangle
    vec3f ss;
    if (mesh->s) {
      ss = (b0 * mesh->s[v[0]] + b1 * mesh->s[v[1]] + b2 * mesh->s[v[2]]);
      if (ss.squared_length() > 0) {
        ss = unit_vector(ss);
      }  else {
        ss = unit_vector(rec.dpdu);
      }
    } else {
      ss = unit_vector(rec.dpdu);
    }
    
    // Compute shading bitangent _ts_ for Triangle and adjust _ss_
    vec3f ts = cross(ss, ns.convert_to_vec3());
    if (ts.squared_length() > 0.f) {
      ts = unit_vector(ts);
      ss = cross(ts, ns.convert_to_vec3());
    } else {
      CoordinateSystem(ns.convert_to_vec3(), &ss, &ts);
    }
    
    // Compute $\dndu$ and $\dndv$ for Triangle shading geometry
    normal3f dndu, dndv;
    if (mesh->n) {
      // Compute deltas for Triangle partial derivatives of normal
      vec2f duv02 = uv[0] - uv[2];
      vec2f duv12 = uv[1] - uv[2];
      normal3f dn1 = mesh->n[v[0]] - mesh->n[v[2]];
      normal3f dn2 = mesh->n[v[1]] - mesh->n[v[2]];
      Float determinant = DifferenceOfProducts(duv02[0], duv12[1], duv02[1], duv12[0]);
      bool degenerateUV = std::abs(determinant) < 1e-8;
      if (degenerateUV) {
        // We can still compute dndu and dndv, with respect to the
        // same arbitrary coordinate system we use to compute dpdu
        // and dpdv when this happens. It's important to do this
        // (rather than giving up) so that ray differentials for
        // rays reflected from Triangles with degenerate
        // parameterizations are still reasonable.
        vec3f dn = cross((mesh->n[v[2]] - mesh->n[v[0]]).convert_to_vec3(),
                         (mesh->n[v[1]] - mesh->n[v[0]]).convert_to_vec3());
        if (dn.squared_length() == 0) {
          dndu = dndv = normal3f(0, 0, 0);
        } else {
          vec3f dnu, dnv;
          CoordinateSystem(dn, &dnu, &dnv);
          dndu = normal3f(dnu);
          dndv = normal3f(dnv);
        }
      } else {
        Float invDet = 1 / determinant;
        dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invDet;
        dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * invDet;
      }
    } else {
      dndu = dndv = normal3f(0, 0, 0);
    }
    if (reverseOrientation) {
      ts = -ts;
    } 
    rec.dpdu = ss;
    rec.dpdv = ts;
    rec.dndu = dndu;
    rec.dndv = dndv;
    
    // isect->SetShadingGeometry(ss, ts, dndu, dndv, true);
  }
  if(bump_tex) {
    point3f bvbu = bump_tex->value(u1, v1, rec.p);
    Float displace = bump_tex->raw_value(u1, v1, rec.p);
    rec.bump_normal = cross(rec.dpdu + bvbu.x() * rec.normal.convert_to_vec3() + displace * rec.dndu.convert_to_vec3(),
                            rec.dpdv - bvbu.y() * rec.normal.convert_to_vec3() + displace * rec.dndv.convert_to_vec3());
      
    rec.has_bump = true;
  } else {
    rec.normal = unit_vector(cross(rec.dpdu, rec.dpdv));
    rec.has_bump = false;
  }
  
  return true;
}


bool Triangle::bounding_box(Float t0, Float t1, aabb& box) const {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  point3f min_v(fmin(fmin(a.x(), b.x()), c.x()), 
                fmin(fmin(a.y(), b.y()), c.y()), 
                fmin(fmin(a.z(), b.z()), c.z()));
  point3f max_v(fmax(fmax(a.x(), b.x()), c.x()), 
                fmax(fmax(a.y(), b.y()), c.y()), 
                fmax(fmax(a.z(), b.z()), c.z()));
  
  point3f difference = max_v + -min_v;
  
  if (difference.x() < 1E-5) max_v.e[0] += 1E-5;
  if (difference.y() < 1E-5) max_v.e[1] += 1E-5;
  if (difference.z() < 1E-5) max_v.e[2] += 1E-5;
  
  box = aabb(min_v, max_v);
  return(true);
}

Float Triangle::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) { 
  hit_record rec;
  if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, rng)) {
    Float distance = rec.t * rec.t * v.squared_length();;
    Float cosine = dot(v, rec.normal);
    return(distance / (cosine * Area()));
  }
  return 0; 
}

Float Triangle::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) { 
  hit_record rec;
  if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, sampler)) {
    Float distance = rec.t * rec.t * v.squared_length();;
    Float cosine = dot(v, rec.normal);
    return(distance / (cosine * Area()));
  }
  return 0; 
}

vec3f Triangle::random(const point3f& origin, random_gen& rng, Float time) {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  Float r1 = rng.unif_rand();
  Float r2 = rng.unif_rand();
  Float sr1 = sqrt(r1);
  point3f random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin); 
}
vec3f Triangle::random(const point3f& origin, Sampler* sampler, Float time) {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  vec2f u = sampler->Get2D();
  Float r1 = u.x();
  Float r2 = u.y();
  Float sr1 = sqrt(r1);
  point3f random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin); 
}

Float Triangle::Area() const {
  // Get triangle vertices in _p0_, _p1_, and _p2_
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  return 0.5 * cross(p1 - p0, p2 - p0).length();
}
