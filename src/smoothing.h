#ifndef SMOOTHINGH
#define SMOOTHINGH

#include "tinyobj/tiny_obj_loader.h"
#include "vec3.h"
#include "mathinline.h"

static void CalcNormal(vec3f& N, vec3f& v0, vec3f& v1, vec3f& v2) {
  vec3f v10 = v1 - v0;
  vec3f v20 = v2 - v0;
  N = cross(v10,v20);
  N.make_unit_vector();
}

// Check if `mesh_t` contains smoothing group id.
bool hasSmoothingGroup(const tinyobj::shape_t& shape) {
  for (size_t i = 0; i < shape.mesh.smoothing_group_ids.size(); i++) {
    if (shape.mesh.smoothing_group_ids[i] > 0) {
      return true;
    }
  }
  return false;
}

void computeSmoothingNormals(const tinyobj::attrib_t& attrib, const tinyobj::shape_t& shape,
                             std::map<int, vec3>& smoothVertexNormals) {
  smoothVertexNormals.clear();
  std::map<int, vec3>::iterator iter;
  
  for (size_t f = 0; f < shape.mesh.indices.size() / 3; f++) {
    // Get the three indexes of the face (all faces are triangular)
    tinyobj::index_t idx0 = shape.mesh.indices[3 * f + 0];
    tinyobj::index_t idx1 = shape.mesh.indices[3 * f + 1];
    tinyobj::index_t idx2 = shape.mesh.indices[3 * f + 2];
    
    // Get the three vertex indexes and coordinates
    int vi[3];      // indexes
    float v[3][3];  // coordinates
    
    for (int k = 0; k < 3; k++) {
      vi[0] = idx0.vertex_index;
      vi[1] = idx1.vertex_index;
      vi[2] = idx2.vertex_index;
      assert(vi[0] >= 0);
      assert(vi[1] >= 0);
      assert(vi[2] >= 0);
      
      v[0][k] = attrib.vertices[3 * vi[0] + k];
      v[1][k] = attrib.vertices[3 * vi[1] + k];
      v[2][k] = attrib.vertices[3 * vi[2] + k];
    }
    
    // Compute the normal of the face
    float normal[3];
    CalcNormal(normal, v[0], v[1], v[2]);
    
    // Add the normal to the three vertexes
    for (size_t i = 0; i < 3; ++i) {
      iter = smoothVertexNormals.find(vi[i]);
      if (iter != smoothVertexNormals.end()) {
        // add
        iter->second.v[0] += normal[0];
        iter->second.v[1] += normal[1];
        iter->second.v[2] += normal[2];
      } else {
        smoothVertexNormals[vi[i]].v[0] = normal[0];
        smoothVertexNormals[vi[i]].v[1] = normal[1];
        smoothVertexNormals[vi[i]].v[2] = normal[2];
      }
    }
    
  }  // f
  
  // Normalize the normals, that is, make them unit vectors
  for (iter = smoothVertexNormals.begin(); iter != smoothVertexNormals.end();
  iter++) {
    make.(iter->second);
  }
  
}  // computeSmoothingNormals
  
#endif
