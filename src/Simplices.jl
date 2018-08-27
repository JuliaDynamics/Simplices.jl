__precompile__(true)

module Simplices

using Distributions

include("barycentric-coordinates.jl")
include("Binary.jl")
include("Circumsphere.jl")
include("complementary.jl")
include("centroids_radii.jl")
include("embed.jl")
include("geometry.jl")
include("heaviside.jl")
include("intersection-of-boundaries.jl")
include("intersection-of-boundaries-loop.jl")
include("init_functions.jl")
include("NullSpace.jl")
include("polytope-generating-vertices.jl")
include("QR.jl")
include("SimplexChecks.jl")
include("shared-vertices.jl")
include("sharing-a-face.jl")
include("some-vertex-in-circumsphere.jl")
include("simplexoperations.jl")
include("simplexintersect.jl")
include("simplex_sampling.jl")
include("simplex_split.jl")
include("simplex_subdivision.jl")
include("simplex_subdivision_single.jl")
include("simplexvolumes.jl")
include("TriangulationPolytopeFaces.jl")
include("TriangulationNonSimplicialFaces.jl")
include("Update.jl")
include("volume-computation.jl")

export BarycentricCoordinates,
        centroid,
        Circumsphere,
        childsimplex,
        childpoint,
        childpoints,
        Delaunay,
        even_sampling, even_sampling_rules, evenly_sample,
        heaviside,
        heaviside0,
        intersecting_simplices,
        intersectingvertices,
        IntersectionOfBoundaries,
        IntersectionOfBoundaries_NoStorage,
        issingular,
        nontrivially_intersecting_simplices,
        nontrivial_intersection,
        outsidepoint,
        orientation,
        PolytopeGeneratingVertices,
        radius,
        SharedVertices,
        SharedFaceVolume,
        SharedFaceVertices,
        shared_vertex_intersection,
        simplices_sharing_vertices,
        simplexintersection,
        SomeVertexInCircumsphere,
        subsample_coeffs,
        tensordecomp,
        volume,
        VolumeComputation

end
