import HDF5: h5open

function savesketch(filename::String, sketch::AbstractSketch)
    h5open(filename, "w") do h5f
        dset = h5f["sketch", "blosc", 9] = getarray(sketch)
        merge!(attrs(dset), h5attrs(sketch))
    end
end

function loadsketch(filename::String)
    h5open(filename, "r") do h5f
        dset = h5f["sketch"]
        attr = attrs(dset)
        sketchtype = attr["sketchtype"]
        sketch = nothing
        if sketchtype == CountMinSketch
            sketch = CountMinSketch{eltype(dset)}(size(dset)...)
        else
            error("Invalid sketch file")
        end
        getarray(sketch) = read(sketch)
        return sketch
    end
end

#function readcms!(filename::AbstractString, cms::CountMinSketch)
#    h5open(filename, "r") do h5f
#        sketch = h5f["sketch"]
#        if eltpye(sketch) != eltype(cms)
#            error("Type mismatch reading sketch: use readcms()")
#        end
#        # Sketch is stored column major
#        ts, nt = size(sketch)
#        cms.sketch = read(sketch)
#        cms.tables = nt
#        cms.tablesize = ts
#        return cms
#    end
#end
