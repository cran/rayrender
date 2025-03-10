% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csg_construct.R
\name{csg_plane}
\alias{csg_plane}
\title{CSG Plane}
\usage{
csg_plane(x = 0, y = 0, z = 0, normal = c(0, 1, 0), width_x = 4, width_z = 4)
}
\arguments{
\item{x}{Default `0`. An x-coordinate on the plane.}

\item{y}{Default `0`. A y-coordinate on the plane.}

\item{z}{Default `0`. A z-coordinate on the plane.}

\item{normal}{Default `c(0,1,0)`. Surface normal of the plane.}

\item{width_x}{Default `10`.}

\item{width_z}{Default `10`.}
}
\value{
List describing the plane in the scene.
}
\description{
Note: This shape isn't closed, so there may be odd lighting issues if it's oriented the wrong
way.
}
\examples{
if(run_documentation()) {
#Generate a plane
csg_object(csg_plane(width_x=4, width_z=4), material=diffuse(checkercolor="purple")) \%>\% 
  add_object(sphere(y=5,x=5,material=light(intensity=40))) \%>\% 
  render_scene(clamp_value=10, samples=16)
 }
if(run_documentation()) {
#Combine the plane with a sphere
csg_object(csg_combine(
    csg_sphere(radius=0.5),
    csg_plane(width_x=4, width_z=4,y=-0.5), 
    operation="blend"),material=diffuse(checkercolor="purple")) \%>\% 
  add_object(sphere(y=5,x=5,material=light(intensity=40))) \%>\% 
  render_scene(clamp_value=10, samples=16)
  }
if(run_documentation()) {
#Re-orient the plane using the normal and 
csg_object(csg_combine(
    csg_sphere(radius=0.5),
    csg_plane(normal = c(1,1,0),width_x=4, width_z=4,y=-0.5), 
    operation="blend"),material=diffuse(checkercolor="purple")) \%>\% 
  add_object(sphere(y=5,x=5,material=light(intensity=40))) \%>\% 
  render_scene(clamp_value=10, samples=16)
}
}
