% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{raymesh_model}
\alias{raymesh_model}
\title{`raymesh` model}
\usage{
raymesh_model(
  mesh,
  x = 0,
  y = 0,
  z = 0,
  flip_transmittance = TRUE,
  verbose = FALSE,
  importance_sample_lights = FALSE,
  calculate_consistent_normals = TRUE,
  override_material = TRUE,
  material = diffuse(),
  angle = c(0, 0, 0),
  order_rotation = c(1, 2, 3),
  flipped = FALSE,
  scale = c(1, 1, 1),
  validate_mesh = TRUE
)
}
\arguments{
\item{mesh}{A `raymesh` object. Pulls the vertex, index, texture coordinates, 
normals, and material information.}

\item{x}{Default `0`. x-coordinate to offset the model.}

\item{y}{Default `0`. y-coordinate to offset the model.}

\item{z}{Default `0`. z-coordinate to offset the model.}

\item{flip_transmittance}{Default `TRUE`. Flips `(1-t)` the transmittance values to match the way the colors
would be interpreted in a rasterizer (where it specifies the transmitted color). Turn off to specify
the attenuation values directly.}

\item{verbose}{Default `FALSE`. If `TRUE`, prints information about the mesh to the console.}

\item{importance_sample_lights}{Default `TRUE`. Whether to importance sample lights specified in the OBJ material
(objects with a non-zero Ke MTL material).}

\item{calculate_consistent_normals}{Default `TRUE`. Whether to calculate consistent vertex normals to prevent energy 
loss at edges.}

\item{override_material}{Default `TRUE`. If `TRUE`, overrides the material specified in the 
`raymesh` object with the one specified in `material`.}

\item{material}{Default  \code{\link{diffuse}}, but ignored unless `override_material = TRUE`. The material, called from one of the material 
functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.}

\item{angle}{Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.}

\item{order_rotation}{Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".}

\item{flipped}{Default `FALSE`. Whether to flip the normals.}

\item{scale}{Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
number, the object will be scaled uniformly.
Note: emissive objects may not currently function correctly when scaled.}

\item{validate_mesh}{Default `TRUE`. Validates the `raymesh` object using `rayvertex::validate_mesh()` 
before parsing to ensure correct parsing. Set to `FALSE` to speed up scene construction if `raymesh_model()` 
is taking a long time (Note: this does not affect rendering time).}
}
\value{
Single row of a tibble describing the raymesh model in the scene.
}
\description{
Load an `raymesh` object, as specified in the `rayvertex` package.
}
\examples{
#Render a simple raymesh object
library(rayvertex)
if(run_documentation()) {
raymesh_model(sphere_mesh(position = c(-1, 0, 0),
              material = material_list(transmittance = "red"))) \%>\%
  add_object(generate_ground(material = diffuse(checkercolor="grey20"))) \%>\%
  render_scene(fov = 30, samples=128, sample_method="sobol_blue")
}

# We create a complex rayvertex mesh, using the `rayvertex::add_shape` function which
# creates a new `raymesh` object out of individual `raymesh` objects
rm_scene = sphere_mesh(position = c(-1, 0, 0),
            material = material_list(transmittance = "red")) \%>\% 
    add_shape(sphere_mesh(position = c(1, 0, 0),
            material = material_list(transmittance = "green", ior = 1.5)))

# Pass the single raymesh object to `raymesh_model()`
# `raymesh_model()`
if(run_documentation()) {
raymesh_model(rm_scene) \%>\%
  add_object(generate_ground(material = diffuse(checkercolor="grey20"))) \%>\%
  render_scene(fov = 30, samples=128, sample_method="sobol_blue")
}

# Set `flip_transmittance = FALSE` argument to specify attenuation coefficients directly
# (as specified in the `dielectric()` material). We change the material's numerical attenuation
# constants using `rayvertex::change_material`
rm_scene_new= change_material(rm_scene, transmittance = c(1,2,0.3), id = 1) \%>\% 
  change_material(transmittance = c(3,1,2), id = 2)
if(run_documentation()) {
raymesh_model(rm_scene_new, flip_transmittance = FALSE) \%>\%
  add_object(generate_ground(material = diffuse(checkercolor="grey20"))) \%>\%
  render_scene(fov = 30, samples=128, sample_method="sobol_blue")
}

# Override the material specified in the `raymesh` object and render the scene
if(run_documentation()) {
raymesh_model(rm_scene,
              material = dielectric(attenuation = "dodgerblue2", attenuation_intensity = 4), 
  override_material = TRUE) \%>\%
  add_object(generate_ground(material = diffuse(checkercolor="grey20"))) \%>\%
  render_scene(fov = 30, samples=128, sample_method="sobol_blue")
}

# Adjusting the scale, position, and rotation parameters of the `raymesh` model
if(run_documentation()) {
raymesh_model(rm_scene,
              x = 0, y = 0.5, z = -1, angle = c(0, 0, 20)) \%>\%
  add_object(generate_ground(material = diffuse(checkercolor="grey20"))) \%>\%
  render_scene(fov = 30,lookat=c(0,0.5,0), samples=128, sample_method="sobol_blue")
}
}
