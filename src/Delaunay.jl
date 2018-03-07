__precompile__()
module Delaunay

using PyCall
const scipyspatial = PyNULL()

function __init__()
    copy!(scipyspatial, pyimport_conda("scipy.spatial", "scipy"))
end

function delaunayn(points::Array{Float64})
    py = scipyspatial[:Delaunay](points)
    indices = Array{Int64, 2}(length(py["simplices"]), size(points, 2) + 1)
    pyarray_to_array!(py["simplices"], indices, Int) 
    return indices .+ 1 # Add 1 to account for base difference in indices
end

function pyarray_to_array!(pyobject, arr, T)
    for i = 1:length(pyobject)
        arr[i, :] = get(pyobject, PyVector{T}, i-1) # i-1 because of Python 0 indexing
    end
end

export delaunayn

end

