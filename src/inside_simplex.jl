"""
    childsimplex(parentsimplex::Array{Float64, 2})

 Generates a random simplex which is entirely contained within `parentsimplex`, which is a (dim+1)-by-dim array.
"""
function childsimplex(parentsimplex::Array{Float64, 2}, dim::Int)
    # Convex expansion coefficients of the random simplex
    dim = size(parentsimplex, 2)
    rs = rand(dim + 1, dim + 1)
    normalised_colsums = 1 ./ sum(rs, 2)
    (normalised_colsums .* rs) * parentsimplex
end


function childpoint(parentsimplex::Array{Float64, 2})
    dim = size(parentsimplex, 2)
    # Random linear combination coefficients
    R = rand(1, dim + 1)

    # Normalise the coefficients so that they sum to one. We can then create the new point
    # as a convex linear combination of the vertices of the parent simplex.
    normalised_coeffs = (1 ./ sum(R, 2)) .* R

    normalised_coeffs * parentsimplex
end


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
