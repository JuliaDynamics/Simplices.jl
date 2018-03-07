"""
    `centroid(simplex::Array{Float64, 2})`

Computes the centroid of a simplex given by `(dim+1)`-by-`dim` array, where each row is a
vertex. Returns the centroid as a vertex (a `1`-by-`dim` two-dimensional array).

"""
function centroid(simplex::Array{Float64, 2})
    n = size(simplex, 1) # Dimension of the space the simplex lives in
    (ones(n, 1)/n).' * simplex
end


"""
    `radius(simplex::Array{Float64, 2}, centroid::Array{Float64, 2})``

Compute radius of a simplex (`(dim+1)`-by-`dim` sized `Array{Float64, 2}` given its centroid (`1`-by-`dim` sized `Array{Float64, 2}`).
"""
function radius(simplex::Array{Float64, 2}, centroid::Array{Float64, 2})

    # Express vertices with respect to origin
    dim = size(simplex, 2)
    dists_to_centroid = broadcast(-, centroid, simplex) # subtract centroid from all vertices
    maximum(sqrt.(sum(dists_to_centroid.^2, 2)))
end


"""
    `radius(s::Array{Float64, 2})`

Compute radius of a simplex `s`, represented by a `Array{Float64, 2}` of size `(dim+1)`-by-`dim`.
"""
function radius(simplex::Array{Float64, 2})

    # Express vertices with respect to origin
    dim = size(simplex, 2)
    dists_to_centroid = broadcast(-, centroid(simplex), simplex) # subtract centroid from all vertices
    maximum(sqrt.(sum(dists_to_centroid.^2, 2)))
end


"""
    orientation(simplex::Array{Float64, 2})

Compute orientation of a `simplex`, represented by a `Array{Float64, 2}` of size `(dim+1)`-by-`dim`.
"""
function orientation(simplex::Array{Float64, 2})
    dim = size(simplex, 2)
    hcat(ones(dim + 1, 1), simplex) |> det
end

"""
    volume(simplex::Array{Float64, 2})

Compute the volume of a `simplex`, represented by a `Array{Float64, 2}` of size `(dim+1)`-by-`dim`.
"""
function volume(simplex::Array{Float64, 2})
    orientation(simplex) |> abs
end

"""
    childsimplex(parentsimplex::Array{Float64, 2})

 Generates a random simplex which is entirely contained within `parentsimplex`, which is a (dim+1)-by-dim array.
"""
function childsimplex(parentsimplex::Array{Float64, 2})
    # Convex expansion coefficients of the random simplex
    dim = size(parentsimplex, 2)
    rs = rand(dim + 1, dim + 1)
    normalised_colsums = 1 ./ sum(rs, 2)
    (normalised_colsums .* rs) * parentsimplex
end

"""
    issingular(simplex::Array{Float64, 2})

Determines if a simplex is singular by checking if any of its vertices are identical.
"""
function issingular(simplex::Array{<:Number, 2})
    size(unique(simplex, 1), 1) != size(simplex, 1)
end

"""
    insidepoints(npts::Int, parentsimplex::Array{Float64, 2})

Generates `npts` points that located inside `parentsimplex`.
"""
function insidepoints(npts::Int, parentsimplex::Array{Float64, 2})
    dim = size(parentsimplex, 2)
    # Random linear combination coefficients
    R = rand(npts, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, 2)) .* R
    normalised_coeffs * parentsimplex
end

"""
    `outsidepoint(parentsimplex::Array{Float64, 2})`

Generate a point that is guaranteed to lie outside `parentsimplex`.
"""
function outsidepoint(parentsimplex::Array{Float64, 2})
    dim = size(parentsimplex, 2)
    # Random linear combination coefficients
    R = rand(1, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, 2)) .* R
    normalised_coeffs[1] += 1

    normalised_coeffs * parentsimplex
end

"""
    outsidepoints(npts::Int, parentsimplex::Array{Float64, 2})

Generates `npts` points that located outside `parentsimplex`.
"""
function outsidepoints(npts::Int, parentsimplex::Array{Float64, 2})
    vcat([outsidepoint(parentsimplex) for i in 1:npts]...)
end



"""
    `nontrivially_intersecting_simplices(dim::Int)`

Genereate a set of non-trivially intersecting `dim`-dimensional simplices (i.e. they don't intersect along boundaries or vertices only).
"""
function nontrivially_intersecting_simplices(dim::Int)
    rs = rand(dim + 1, dim) # random simplex with dim+1 vertices
    n_inside = rand(1:dim, 1)[1]

    verts_inside = insidepoints(n_inside, rs)

    verts_outside = vcat([outsidepoint(rs) for i in 1:(dim+1-n_inside)]...)

    return rs, vcat(verts_inside, verts_outside)
end


"""
    `nontrivially_intersecting_simplices(dim::Int)`

Genereate a set of non-trivially intersecting `dim`-dimensional simplices (i.e. they don't intersect along boundaries or vertices only).
"""
function simplices_sharing_vertices(dim::Int)
    rs = rand(dim + 1, dim) # random simplex with dim+1 vertices

    # Determine how many vertices should be shared and how many are remaining
    n_shared_verts = rand(1:dim, 1)[1]
    n_remaining_verts = dim + 1 - n_shared_verts

    # Pick the shared vertices and generate remaining vertices
    new_simplex = zeros(Float64, 0, dim)
    shared_verts = rs[rand(1:n_shared_verts, n_shared_verts), :]

    if ndims(shared_verts) == 1
        new_simplex = vcat(new_simplex, shared_verts.')
    else
        new_simplex = vcat(new_simplex, shared_verts)
    end
    for i = 1:n_remaining_verts
        if rand(Bool) == true
            new_simplex  = vcat(new_simplex, insidepoints(1, rs))
        else
            new_simplex  = vcat(new_simplex, outsidepoints(1, rs))
        end
    end

    return rs, new_simplex
end

"""
    intersecting_simplices(;dim::Int = 3, intersection_type = "nontrivial")

Generate a pair of `dim` dimensional intersecting simplices. Each resulting
simplex is a `(dim+1)`-by-`dim` two-dimensional array. The `intersection_type` argument can be either `"nontrivial"` or `"sharingvertices"`.
"""
function intersecting_simplices(;dim::Int = 3, intersection_type = "nontrivial")
    if intersection_type == "nontrivial"
        return nontrivially_intersecting_simplices(dim)
    elseif intersection_type == "sharingvertices"
        return simplices_sharing_vertices(dim)
    end
end
