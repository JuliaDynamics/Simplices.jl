"""
    childpoint(parentsimplex::AbstractArray{Float64, 2})

Generates a random simplex which is entirely contained within `parentsimplex`,
which is a (dim+1)-by-dim array.
"""
function childpoint(parentsimplex::AbstractArray{Float64, 2})
    dim = size(parentsimplex, 2)
    # Random linear combination coefficients
    R = rand(Float64, 1, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create
	# the new point as a convex linear combination of the vertices of the parent
	# simplex.
    R = (1 ./ sum(R, 2)) .* R

    R * parentsimplex
end

"""
    childpoint(parentsimplex::AbstractArray{Float64, 2})

Generates a random simplex which is entirely contained within `parentsimplex`,
which is a (dim+1)-by-dim array.
"""
function childpoint(parentsimplex::AbstractArray{Float64, 2}, dim::Int)
    # Random linear combination coefficients
    R = rand(Float64, 1, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create
	# the new point as a convex linear combination of the vertices of the parent
	# simplex.
    R = (1 ./ sum(R, 2)) .* R

    R * parentsimplex
end


"""
Create `n` (default: n = 100) points that lies inside `parent_simplex`,
which is a (dim+1)-by-dim array.
"""
function childpoints(parent_simplex::AbstractArray{Float64, 2}, n::Int)
    dim = size(parent_simplex, 2)

    # Random linear combination coefficients
    R = rand(Float64, n, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create
    # the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, 2)) .* R

    normalised_coeffs * parent_simplex
end


"""
Create a point that lies inside `parent_simplex`, which is a (dim+1)-by-dim
array.
"""
function outsidepoint(parentsimplex::AbstractArray{Float64, 2})
    dim = size(parentsimplex, 2)
    # Random linear combination coefficients
    R = rand(Float64, 1, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, 2)) .* R
    normalised_coeffs[1] += 1

    normalised_coeffs * parentsimplex
end
