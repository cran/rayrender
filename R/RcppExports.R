# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

PrintClassSizes <- function() {
    invisible(.Call(`_rayrender_PrintClassSizes`))
}

has_gui_capability <- function() {
    .Call(`_rayrender_has_gui_capability`)
}

cppdef_HAS_OIDN <- function() {
    .Call(`_rayrender_cppdef_HAS_OIDN`)
}

cppdef_HAS_NEON <- function() {
    .Call(`_rayrender_cppdef_HAS_NEON`)
}

cppdef_HAS_SSE <- function() {
    .Call(`_rayrender_cppdef_HAS_SSE`)
}

cppdef_HAS_SSE2 <- function() {
    .Call(`_rayrender_cppdef_HAS_SSE2`)
}

cppdef_HAS_SSE3 <- function() {
    .Call(`_rayrender_cppdef_HAS_SSE3`)
}

cppdef_HAS_SSE41 <- function() {
    .Call(`_rayrender_cppdef_HAS_SSE41`)
}

render_animation_rcpp <- function(scene, camera_info, scene_info, render_info, camera_movement, start_frame, end_frame, filenames, post_process_frame, toneval, bloom, write_image, transparent_background) {
    invisible(.Call(`_rayrender_render_animation_rcpp`, scene, camera_info, scene_info, render_info, camera_movement, start_frame, end_frame, filenames, post_process_frame, toneval, bloom, write_image, transparent_background))
}

render_scene_rcpp <- function(scene, camera_info, scene_info, render_info) {
    .Call(`_rayrender_render_scene_rcpp`, scene, camera_info, scene_info, render_info)
}

tonemap_image <- function(routput, goutput, boutput, toneval) {
    .Call(`_rayrender_tonemap_image`, routput, goutput, boutput, toneval)
}

