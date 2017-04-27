"""
A Count-min Sketch
"""
type CountMinSketch{T<:Unsigned} <: AbstractSketch
    # We store the sketch as a Matrix
    sketch::Matrix{T}
    seeds::Vector{UInt}

    """
    The constructor for `CountMinSketch`es

    Arguments
    ---------

    * `tables::Int`: The number of tables to create
    * `tablesize::Int`: The size of each table
    """
    function CountMinSketch(array::Matrix{T})
        tables, tablesize = size(array)
        if !(1 <= tables <= 20)
            error("Must have between 1 and 20 tables")
        end
        tablesize > 1 || error("Table size must be greater than 1")
        seeds = map(hash, 1:tables)
        return new(array, seeds)
    end
    function CountMinSketch(tables::Int, tablesize::Int)
        if !(1 <= tables <= 20)
            error("Must have between 1 and 20 tables")
        end
        tablesize > 1 || error("Table size must be greater than 1")
        sketch = zeros(T, tables, tablesize)
        seeds = map(hash, 1:tables)
        return new(sketch, seeds)
    end
end



"""
Adds `count` `item`s to a CountMinSketch

This function will never allow the counts for an item to overflow or underflow.

Arguments
---------
* `cms`: A count-min sketch
* `item`: any hashable item
* `count::Integer`: a (potentially negative) number of `item`s to add.
"""
function add!(cms::CountMinSketch, item, count::Int)
    for i in 1:size(cms)[1]
        offset = (hash(item, cms.seeds[i]) % size(cms)[2]) + 1
        @inbounds x = cms.sketch[i, offset]
        @inbounds cms.sketch[i, offset] = clamp(x + count,
                                                typemin(eltype(cms)),
                                                typemax(eltype(cms)))
    end
end


"""
Determine count of an item in a CountMinSketch.

NB: This function returns an estimate of the true count. This estimate is never
    lower than the true count, but has a very low probablity of being higher.

Arguments
---------
* `cms`: A count-min sketch
* `item`: any hashable item
"""
function Base.getindex{T}(cms::CountMinSketch{T}, item)
    minval = typemax(T)
    for i in 1:size(cms)[1]
        offset = (hash(item, cms.seeds[i]) % size(cms)[2]) + 1
        @inbounds val = cms.sketch[i, offset]
        minval = val < minval ? val : minval
    end
    return minval
end



"""
Element type of a CountMinSketch
"""
function Base.eltype(cms::CountMinSketch)
    return eltype(cms.sketch)
end


"""
Shape of a CountMinSketch, i.e. (number of tables, tablesize).
"""
function Base.size(cms::CountMinSketch)
    return size(cms.sketch)
end


"""
Estimates the collision or aliasing rate of a `CountMinSketch`.

The collision rate is the probability that a random key will be added in such a
way that it is indistinguishable to another random key. This is equivalent to
the product of the probablities of colliding in each table independently, which
is simply the occupancy rate of each table.

Generally avoid using a `CountMinSketch` with a collision rate above about 0.05
- 0.10.
"""
function collisionrate(cms::CountMinSketch)
    occupancy = sum(cms.sketch .> 0, 2)
    # The total FPR is the product of all rates, as we assume each is truly
    # independent
    rate = prod(float(occupancy) / float(size(cms)[2]))
    return rate
end
