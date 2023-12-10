// #include <Rcpp.h>
// #include <memory>
// 
// #include "float.h"
// #include "vec3.h"
// #include "vec2.h"
// #include "point3.h"
// #include "point2.h"
// #include "normal.h" 
// #include "RayMatrix.h"
// #include "mathinline.h"
// #include "transform.h"
// #include "transformcache.h"
// #include "camera.h"
// #include "float.h"
// #include "buildscene.h"
// #include "RProgress.h"
// #include "rng.h"
// #include "tonemap.h"
// #include "infinite_area_light.h"
// #include "adaptivesampler.h"
// #include "sampler.h"
// #include "color.h"
// #include "integrator.h"
// #include "debug.h"
// using namespace Rcpp;
// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::depends(RcppThread)]]
// #include "RcppThread.h"
// #include "PreviewDisplay.h"
// 
// using namespace Rcpp;
// 
// class SceneRenderer {
//   // Class members
//   List camera_info;
//   List scene_info;
//   
//   TransformCache transformCache;
//   TransformCache transformCacheBg;
//   
//   RayMatrix routput;
//   RayMatrix goutput;
//   RayMatrix boutput;
//   
//   std::unique_ptr<RayCamera> cam;
//   
//   int nx1, ny1, nn1;
//   
//   std::vector<Float* > textures;
//   std::vector<int* > nx_ny_nn;
//   
//   std::vector<unsigned char * > alpha_textures;
//   std::vector<int* > nx_ny_nn_alpha;
//   
//   std::vector<unsigned char * > bump_textures;
//   std::vector<int* > nx_ny_nn_bump;
//   
//   std::vector<unsigned char * > roughness_textures;
//   std::vector<int* > nx_ny_nn_roughness;
//   
//   std::vector<std::shared_ptr<material> >* shared_materials;
//   
// public:
//   SceneRenderer(List camera_info, List scene_info) 
//     : camera_info(camera_info), scene_info(scene_info), 
//       shared_materials(new std::vector<std::shared_ptr<material> >) {
//     
//   }
//   
//   //Unpack scene info to an Rcpp type
//   void UnpackSceneInfo() {
//     //...
//   }
//   
//   //Unpack Camera Info
//   void UnpackCameraInfo() {
//     //...
//   }
//   
//   //Initialize transformation cache
//   void InitTransformCache() {
//     //...
//   }
//   
//   //Initialize output matrices
//   void InitOutputMatrices() {
//     //...
//   }
//   
//   //Initialize textures
//   void InitTextures() {
//     //...
//   }
//   
//   // Load Textures
//   void LoadTextures() {
//     //...
//   }
//   
//   // Build Scene BVH
//   std::shared_ptr<hitable> BuildSceneBVH() {
//     //...
//   }
//   
//   //Handle Background
//   void HandleBackground() {
//     //...
//   }
//   
//   //Final Steps
//   void FinalSteps() {
//     //...
//   }
// };
// 
// using RendererPtr = std::shared_ptr<SceneRenderer>;
// 
// // [[Rcpp::export]]
// SEXP createRenderer() {
//   RendererPtr renderer(new SceneRenderer());
//   return Rcpp::XPtr<RendererPtr>(new RendererPtr(renderer));
// }
// 
// // [[Rcpp::export]]
// void addElement(SEXP rendererPtr, Element element) {
//   Rcpp::XPtr<RendererPtr> r(rendererPtr);
//   (*r)->addElement(element);
// }
