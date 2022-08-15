#include "trianglemesh.h"

TriangleMesh::TriangleMesh(const Transform &ObjectToWorld,
                           int nTriangles, const int *vertexIndices, int nVertices,
                           const point3f *P, const vec3f *S, const normal3f *N,
                           const point2f *UV,
                           const std::shared_ptr<alpha_texture> &alphaMask)
  : nTriangles(nTriangles), nVertices(nVertices), 
    vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles),
    alphaMask(alphaMask) {
  p.reset(new point3f[nVertices]);
  for (int i = 0; i < nVertices; ++i) {
    p[i] = ObjectToWorld(P[i]);
  }
  
  if (UV) {
    uv.reset(new point2f[nVertices]);
    memcpy(uv.get(), UV, nVertices * sizeof(point2f));
  }
  if (N) {
    n.reset(new normal3f[nVertices]);
    for (int i = 0; i < nVertices; ++i) {
      n[i] = ObjectToWorld(N[i]);
    }
  }
  if (S) {
    s.reset(new vec3f[nVertices]);
    for (int i = 0; i < nVertices; ++i) {
      s[i] = ObjectToWorld(S[i]);
    }
  }
}
