#ifndef INTERACTIONH;
#define INTERACTIONH;

struct Interaction {
  Interaction() : time(0) { }
  Interaction(const vec3 &p, const vec3 &n, const vec3 &pError,
              const Vector3f &wo, Float time,
              const MediumInterface &mediumInterface)
    : p(p), time(time),  wo(wo), n(n),
      mediumInterface(mediumInterface) { }
  bool IsSurfaceInteraction() const {
    return n != Normal3f();
  }
  Ray SpawnRay(const Vector3f &d) const {
    Point3f o = OffsetRayOrigin(p, n, d);
    return Ray(o, d, Infinity, time, GetMedium(d));
  }
  Ray SpawnRayTo(const Point3f &p2) const {
    Point3f origin = OffsetRayOrigin(p, pError, n, p2 - p);
    Vector3f d = p2 - origin;
    return Ray(origin, d, 1 - ShadowEpsilon, time, GetMedium(d));
  }
  Ray SpawnRayTo(const Interaction &it) const {
    Point3f origin = OffsetRayOrigin(p, pError, n, it.p - p);
    Point3f target = OffsetRayOrigin(it.p, it.pError, it.n, origin - it.p);
    Vector3f d = target - origin;
    return Ray(origin, d, 1-ShadowEpsilon, time, GetMedium(d));
  }
  Interaction(const Point3f &p, const Vector3f &wo, Float time,
              const MediumInterface &mediumInterface)
    : p(p), time(time), wo(wo), mediumInterface(mediumInterface) { }
  Interaction(const Point3f &p, Float time,
              const MediumInterface &mediumInterface)
    : p(p), time(time), mediumInterface(mediumInterface) { }
  bool IsMediumInteraction() const { return !IsSurfaceInteraction(); }
  const Medium *GetMedium(const Vector3f &w) const {
    return Dot(w, n) > 0 ? mediumInterface.outside :
    mediumInterface.inside;
  }
  const Medium *GetMedium() const {
    Assert(mediumInterface.inside == mediumInterface.outside);
    return mediumInterface.inside;
  }
  
  <<Interaction Public Data>> 
  Point3f p;
  Float time;
  Vector3f pError;
  Vector3f wo;
  Normal3f n;
  MediumInterface mediumInterface;
    
};

#endif
