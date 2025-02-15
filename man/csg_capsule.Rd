% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csg_construct.R
\name{csg_capsule}
\alias{csg_capsule}
\title{CSG Capsule}
\usage{
csg_capsule(start = c(0, 0, 0), end = c(0, 1, 0), radius = 1)
}
\arguments{
\item{start}{Default `c(0, 0, 0)`. Start point of the capsule, specifying `x`, `y`, `z`.}

\item{end}{Default `c(0, 1, 0)`. End point of the capsule, specifying `x`, `y`, `z`.}

\item{radius}{Default `1`. Capsule radius.}
}
\value{
List describing the capsule in the scene.
}
\description{
CSG Capsule
}
\examples{
if(run_documentation()) {
#Generate a basic capsule:
generate_ground(material=diffuse(checkercolor="grey20")) \%>\% 
  add_object(csg_object(csg_capsule(radius=0.5),material=glossy(color="red"))) \%>\% 
  render_scene(clamp_value=10, samples=16,fov=20)
  }
if(run_documentation()) {
#Change the orientation by specifying a start and end
generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) \%>\% 
  add_object(csg_object(csg_capsule(start = c(-1,0.5,-2), end = c(1,0.5,-2),
  radius=0.5),material=glossy(checkercolor="red"))) \%>\% 
  render_scene(clamp_value=10, samples=16,fov=20,
               lookat=c(0,0.5,-2),lookfrom=c(3,3,10))
 }
if(run_documentation()) {
#Show the effect of changing the radius
generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) \%>\% 
  add_object(csg_object(
    csg_combine(
    csg_capsule(start = c(-1,0.5,-2), end = c(1,0.5,-2), radius=0.5),
    csg_capsule(start = c(-0.5,1.5,-2), end = c(0.5,1.5,-2), radius=0.25)),
    material=glossy(checkercolor="red"))) \%>\% 
  render_scene(clamp_value=10, samples=16,fov=20,
               lookat=c(0,0.5,-2),lookfrom=c(-3,3,10))
   }
if(run_documentation()) {            
#Render a capsule in a Cornell box
generate_cornell() \%>\% 
  add_object(csg_object(
    csg_capsule(start = c(555/2-100,555/2,555/2), end = c(555/2+100,555/2,555/2), radius=100),
    material=glossy(color="dodgerblue4"))) \%>\% 
  render_scene(clamp_value=10, samples=16)
}
}
