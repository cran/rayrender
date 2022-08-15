#ifndef TRIANGLEMESHH
#define TRIANGLEMESHH

#include "texture.h"
#include "transform.h"

struct TriangleMesh {
  TriangleMesh(const Transform &ObjectToWorld, int nTriangles, const int *vertexIndices,
               int nVertices, const point3f *P, const vec3f *S, const normal3f *N,
               const point2f *uv, const std::shared_ptr<alpha_texture> &alphaMask);
  
  const int nTriangles, nVertices;
  std::vector<int> vertexIndices;
  std::unique_ptr<point3f[]> p;
  std::unique_ptr<normal3f[]> n;
  std::unique_ptr<vec3f[]> s;
  std::unique_ptr<point2f[]> uv;
  std::shared_ptr<alpha_texture> alphaMask;
  
};

#endif
