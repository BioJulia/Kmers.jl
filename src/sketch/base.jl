"Base of all sketching data structures"
abstract type AbstractSketch end

"""
Add one `item` to a sketch
"""
function Base.push!(sketch::AbstractSketch, item)
    add!(sketch, item, 1)
end


"""
Remove one `item` from a AbstractSketch
"""
function Base.pop!(sketch::AbstractSketch, item)
    add!(sketch, item, -1)
end

"""
Test for presence of an item in a AbstractSketch.
"""
function Base.haskey(sketch::AbstractSketch, item)
    return sketch[item] > 0
end

function Base.show(io::IO, sketch::AbstractSketch)
    print("$(typeof(sketch)) with shape $(shape(sketch)) and eltype $(eltype(sketch)")
end
