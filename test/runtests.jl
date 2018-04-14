using Simplices
using Base.Test

# installed = Pkg.installed()
# if !("SimplexSplitting" in keys(installed))
#     Pkg.clone("https://github.com/kahaaga/SimplexSplitting.jl")
# else
#     using SimplexSplitting
# end
#
# using PyCall, Conda
# Conda.add("scipy")
# ENV["PYTHON"]= ""; Pkg.build("PyCall"); using PyCall


include("even_sampling.jl")
include("test_delaunay.jl")
include("test_intersectionvol.jl")
include("tensordecomposition.jl")
