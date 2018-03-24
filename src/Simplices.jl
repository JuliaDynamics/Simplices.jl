__precompile__(true)

module Simplices

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
include("inside_simplex.jl")
include("NullSpace.jl")
include("polytope-generating-vertices.jl")
include("QR.jl")
include("SimplexChecks.jl")
include("shared-vertices.jl")
include("sharing-a-face.jl")
include("some-vertex-in-circumsphere.jl")
include("simplexoperations.jl")
include("simplexintersect.jl")
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
        volume,
        VolumeComputation

        # # Simplex splitting stuff
        # centroids_radii2,
        # embed,
        # Embedding,
        # embedding,
        # embedding_ex,
        # embedding_example,
        # example_triangulation,
        # gaussian_embedding,
        # gaussian_embedding_arr,
        # invariantset,
        # refine_triangulation,
        # refine_triangulation_images,
        # refine_recursive,
        # refine_recursive_images,
        # refine_t!,
        # refine_variable_k!,
        # refine_variable_k_new,
        # refine_variable_k_newnew,
        # simplex_split,
        # simplicial_subdivision_single,
        # simplicial_subdivision,
        # simplex_volumes,
        # tensordecomposition,
        # triangulate,
        # Triangulation,
        # triang_from_embedding
        # query_refinement

end
