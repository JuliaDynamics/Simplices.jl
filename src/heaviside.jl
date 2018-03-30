heaviside{T<:Number}(x::T) = ifelse(x <= 0, zero(x), ifelse(x > 0, one(x), oftype(x,0)))
heaviside0{T<:Number}(x::T) = ifelse(x < 0, zero(x), ifelse(x >= 0, one(x), oftype(x, 0)))

func_name1 = Symbol(string("heaviside"))

# generate heaviside functions for variably sized arrays
for t in 1:5
    @eval function $(func_name1){T<:Number}(v::AbstractArray{T, $t})
        h = similar(v)
        for i in eachindex(v)
            h[i] = heaviside(v[i])
        end
        return(h)
    end
end

@eval function $(func_name1){T<:Number}(v::Vector{T})
    h = similar(v)
    for i in eachindex(v)
        h[i] = heaviside(v[i])
    end
    return(h)
end

@eval function $(func_name1){T<:Number}(v::RowVector{T})
    h = similar(v)
    for i in eachindex(v)
        h[i] = heaviside(v[i])
    end
    return(h)
end

func_name2 = Symbol(string("heaviside0"))

# generate heaviside0 functions for variably sized arrays
for t in 1:5
    @eval function $(func_name2){T<:Number}(v::AbstractArray{T, $t})
        h = similar(v)
        for i in eachindex(v)
            h[i] = heaviside0(v[i])
        end
        return(h)
    end
end

@eval function $(func_name2){T<:Number}(v::Vector{T})
    h = similar(v)
    for i in eachindex(v)
        h[i] = heaviside0(v[i])
    end
    return(h)
end

@eval function $(func_name2){T<:Number}(v::RowVector{T})
    h = similar(v)
    for i in eachindex(v)
        h[i] = heaviside0(v[i])
    end
    return(h)
end
