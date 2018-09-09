function complementary(v::Array{Float64}, n::Int64)
  ns = collect(1:n)
  setdiff(1:n, v)
end

function complementary(v::Vector{Float64}, n::Int64)
  ns = collect(1:n)
  setdiff(1:n, v)
end

function complementary(v::Vector{Int}, n::Int64)
  ns = collect(1:n)
  setdiff(1:n, v)
end

function complementary(v::Int64, n::Int64)
  ns = collect(1:n)
  setdiff(1:n, v)
end


function complementary(v::Array{Int64}, n::Int64)
  ns = collect(1:n)
  setdiff(1:n, v)
end
