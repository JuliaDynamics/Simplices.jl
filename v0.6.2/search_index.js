var documenterSearchIndex = {"docs":
[{"location":"simplexintersection/#How-are-intersections-computed?","page":"How to calculate intersections?","title":"How are intersections computed?","text":"","category":"section"},{"location":"simplexintersection/","page":"How to calculate intersections?","title":"How to calculate intersections?","text":"The intersection between two n-dimensional simplices is calculated boundary triangulation. Intersections are computed as follows.","category":"page"},{"location":"simplexintersection/","page":"How to calculate intersections?","title":"How to calculate intersections?","text":"Find minimal set of points generating the intersection volume. These points form a convex polytope Pᵢ.\nTriangulate the faces of Pᵢ into simplices.\nThen combine each resulting boundary simplex with an interior point in Pᵢ. The set of all such combinations now form a triangulation of Pᵢ.\nCalculate the volume of each simplex in the resulting triangulation. Summing over these volumes given the volume of the intersection.","category":"page"},{"location":"funcs/#Function-reference","page":"Function reference","title":"Function reference","text":"","category":"section"},{"location":"funcs/#Intersection-between-simplices","page":"Function reference","title":"Intersection between simplices","text":"","category":"section"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"simplexintersection(S1::Array{Float64, 2}, S2::Array{Float64, 2}; tolerance::Float64 = 1e-10)","category":"page"},{"location":"funcs/#Simplices.simplexintersection-Tuple{Matrix{Float64}, Matrix{Float64}}","page":"Function reference","title":"Simplices.simplexintersection","text":"  simplexintersection(s1::Array{Float64, 2}, s2::Array{Float64, 2};\n    tolerance::Float64 = 1/10^10) -> Float64\n\nComputes the volume of intersection between two n-dimensional simplices by boundary triangulation. The simplices s1 and s2 are arrays of size (dim, dim+1), where each column is a vertex.\n\nNote: the returned volume is not corrected. It should be divided by factorial(dim) to obtain the true volume.\n\nHow are intersections computed?\n\nIntersections are computed as follows:\n\nFind minimal set of points generating the intersection volume. These points form a convex polytope Pᵢ.\nTriangulate the faces of Pᵢ into simplices. \nCombine each boundary simplex with an interior point in Pᵢ. The set of all such combinations form a triangulation of Pᵢ.\nCalculate the volume of each simplex in the resulting triangulation. The sum of these volumes is the volume of the intersection.\n\n\n\n\n\n","category":"method"},{"location":"funcs/#Generating-points-inside/outside-simplex","page":"Function reference","title":"Generating points inside/outside simplex","text":"","category":"section"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"insidepoints(npts::Int, parentsimplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.insidepoints-Union{Tuple{T}, Tuple{Int64, Matrix{T}}} where T<:Number","page":"Function reference","title":"Simplices.insidepoints","text":"insidepoints(npts::Int, parentsimplex::AbstractArray{T, 2}) where {T<:Number}\n\nGenerates npts points that located inside parentsimplex.\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"outsidepoint(parentsimplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.outsidepoint-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Function reference","title":"Simplices.outsidepoint","text":"outsidepoint(parentsimplex::AbstractArray{T, 2}) where {T<:Number}\n\nGenerate a single point that is guaranteed to lie outside parentsimplex.\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"outsidepoints(npts::Int, parentsimplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.outsidepoints-Union{Tuple{T}, Tuple{Int64, Matrix{T}}} where T<:Number","page":"Function reference","title":"Simplices.outsidepoints","text":"outsidepoints(npts::Int, parentsimplex::AbstractArray{T, 2}) where T <: Number\n\nGenerates npts points that located outside parentsimplex.\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"childsimplex(parentsimplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.childsimplex-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Function reference","title":"Simplices.childsimplex","text":"childsimplex(parentsimplex::AbstractArray{T, 2}) where {T<:Number} -> Array{Float64, 2}\n\nGenerates a random simplex which is entirely contained within parentsimplex, which is a (dim+1)-by-dim array.\n\n\n\n\n\n","category":"method"},{"location":"funcs/#Generating-simplices-that-intersect","page":"Function reference","title":"Generating simplices that intersect","text":"","category":"section"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"There are many ways simplices possibly can intersect, but they all boil down to three cases: 1) there is no intersection, 2) they intersect along a common vertex or boundary, or 3) the intersection is more complex. The following functions generate simplices that either share at least one vertex, or simplices that intersect in nontrivial ways (i.e. intersection is not along a vertex or an edge).","category":"page"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"simplices_sharing_vertices(dim::Int)","category":"page"},{"location":"funcs/#Simplices.simplices_sharing_vertices-Tuple{Int64}","page":"Function reference","title":"Simplices.simplices_sharing_vertices","text":"simplices_sharing_vertices(dim::Int) -> Array{Float64, 2}\n\nGenereate a set of non-trivially intersecting dim-dimensional simplices (i.e. they don't intersect along boundaries or vertices only).\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"nontrivially_intersecting_simplices(dim::Int)","category":"page"},{"location":"funcs/#Simplices.nontrivially_intersecting_simplices-Tuple{Int64}","page":"Function reference","title":"Simplices.nontrivially_intersecting_simplices","text":"nontrivially_intersecting_simplices(dim::Int) -> Array{Float64, 2}\n\nGenereate a set of non-trivially intersecting dim-dimensional simplices (i.e. they don't intersect along boundaries or vertices only).\n\n\n\n\n\n","category":"method"},{"location":"funcs/#Properties-of-simplices","page":"Function reference","title":"Properties of simplices","text":"","category":"section"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"orientation(simplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.orientation-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Function reference","title":"Simplices.orientation","text":"orientation(simplex::AbstractArray{T, 2}) where {T<:Number} -> Float64\n\nCompute orientation of a simplex, represented by a Array{Float64, 2} of size (dim+1)-by-dim.\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"volume(simplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.volume-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Function reference","title":"Simplices.volume","text":"volume(simplex::AbstractArray{T, 2}) where {T<:Number} -> Float64\n\nCompute the volume of a simplex, represented by a Array{Float64, 2} of size (dim+1)-by-dim.\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"radius(simplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.radius-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Function reference","title":"Simplices.radius","text":"radius(s::AbstractArray{T, 2}) where {T<:Number} -> Float64\n\nCompute radius of a simplex s, represented by a Array{Float64, 2} of size (dim+1)-by-dim.\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"radius(simplex::Array{T, 2}, centroid::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.radius-Union{Tuple{T}, Tuple{Matrix{T}, Matrix{T}}} where T<:Number","page":"Function reference","title":"Simplices.radius","text":"radius(simplex::Array{T, 2}, centroid::Array{T, 2}) where {T<:Number}\n\nCompute radius of a simplex ((dim+1)-by-dim sized Array{Float64, 2} given its centroid (1-by-dim sized Array{Float64, 2}).\n\n\n\n\n\n","category":"method"},{"location":"funcs/","page":"Function reference","title":"Function reference","text":"centroid(simplex::Array{T, 2}) where {T<:Number}","category":"page"},{"location":"funcs/#Simplices.centroid-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Function reference","title":"Simplices.centroid","text":"centroid(simplex::AbstractArray{Float64, 2}) where {T<:Number} -> Array{Float64, 2}\n\nComputes the centroid of a simplex given by (dim+1)-by-dim array, where each row is a vertex. Returns the centroid as a vertex (a 1-by-dim two-dimensional array).\n\n\n\n\n\n","category":"method"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"We'll use a few of the functions used for testing the package to demonstrate its usage. These are","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"simplices_sharing_vertices. Generate a set of simplices which intersect in some arbitrary way, but sharing at least one vertex.\nnontrivially_intersecting_simplices. Generate a set of non-trivially intersecting simplices (i.e. intersections are not only along boundaries or vertices).\nchildsimplex. Generate a simplex completely contained within a parent simplex.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Note that these functions take as inputs simplices of shape (dim + 1, dim). This will be fixed in a future release.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In the examples, some of the functions used to generate simplices return the simplex arrays with rows as columns. Therefore, we transpose the simplices before calling simplexintersection.","category":"page"},{"location":"examples/#Simplices-sharing-vertices","page":"Examples","title":"Simplices sharing vertices","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Let's compute the intersection between a set of simplices that share at least one vertex.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Simplices\ns₁, s₂ = transpose.(simplices_sharing_vertices(3))\n\nsimplexintersection(s₁, s₂)","category":"page"},{"location":"examples/#Nontrivially-intersecting-simplices","page":"Examples","title":"Nontrivially intersecting simplices","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Simplices can also intersect in nontrivial ways, meaning that they have an  intersection beyond a common boundary or vertex.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Simplices\ns₁, s₂ = transpose.(nontrivially_intersecting_simplices(3))\nsimplexintersection(s₁, s₂)","category":"page"},{"location":"examples/#One-simplex-fully-contained-within-the-other","page":"Examples","title":"One simplex fully contained within the other","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"We'll generate a random simplex s₁, then generate a simplex s₂ fully contained within that simplex. If s₂ is fully contained, the intersection volume should be the volume of s₂.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Simplices\n\nDs =  2:7;\nintersection_vols = zeros(Float64, length(Ds));\nanalytical_vols   = zeros(Float64, length(Ds));\nfor i = 1:length(Ds)\n    s₁ = rand(Ds[i] + 1, Ds[i])\n    s₂ = childsimplex(s₁)\n    intersection_vols[i] = simplexintersection(transpose(s₁), transpose(s₂))\n    analytical_vols[i] = volume(s₂)\nend\n\n# Within numerical error, the results should be the same.\nall([isapprox(intersection_vols[i], analytical_vols[i]; atol = 1e-9) \n    for i = 1:length(Ds)])","category":"page"},{"location":"examples/#Simplices-are-identical","page":"Examples","title":"Simplices are identical","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"If simplices are identical, the intersection volume should equal the volume of either simplex:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Simplices\n\ns₁ = rand(4, 3); s₂ = s₁;\n\nsimplexintersection(transpose(s₁), transpose(s₂)) .≈ volume(s₁) .≈ volume(s₂)","category":"page"},{"location":"#Simplices.jl-documentation","page":"Overview","title":"Simplices.jl documentation","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Julia package for computing exact intersection volumes between n-dimensional simplices.","category":"page"},{"location":"#Usage","page":"Overview","title":"Usage","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Simplex intersections are computed by calling simplexintersection with the simplices in question. The simplices must be arrays of size (dim, dim + 1), so that each vertex of the simplex is a column vector.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Simplex properties, e.g. the radius or centroid, can also be computed for individual simplices. There are also some  functions that can be used to generate pairs of simplices that overlap in certain ways, or points that lie either outside or inside a simplex (insidepoints and outsidepoints).","category":"page"},{"location":"#D-example","page":"Overview","title":"3D example","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"using Simplices\ns₁, s₂ = rand(3, 4), rand(3, 4)\nsimplexintersection(s₁, s₂)","category":"page"},{"location":"#D-example-2","page":"Overview","title":"5D example","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"using Simplices\ns₃, s₄ = rand(4, 5), rand(4, 5)\nsimplexintersection(s₃, s₄)","category":"page"}]
}