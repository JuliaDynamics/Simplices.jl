"""
    Circumsphere(simplex)

Compute the circumsphere of a simplex (a matrix whose columns represent the vertices of the
simplex). Returns a column vector where the first entry is the radius of the circumsphere
and the remaining entries represent the centroid.

"""
function Circumsphere(simplex::Array{Float64, 2})
    n::Int = size(simplex, 1) # Dimension of the space the simplex lives in
    centroid::Array{Float64, 2} = simplex * ones(n + 1, 1) / (n + 1)
    centroidmatrix::Array{Float64, 2} = repmat(centroid, 1, n + 1)

    radius::Float64 = sqrt(maximum(ones(1, n) * ((simplex - centroidmatrix).^2)))

    return [radius; centroid]
end
