using Simplices
using Base.Test


include("test_delaunay.jl")
include("test_intersectionvol.jl")

delaunaytest(10, 3, 10)

speed_test_nontrivial(3, 1)
speed_test_nontrivial(4, 1)
speed_test_nontrivial(5, 1)

speed_test_sharing(3, 1)
speed_test_sharing(4, 1)
speed_test_sharing(5, 1)
