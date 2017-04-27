function intempdir(fn::Function, parent=tempdir())
    tmpdir = mktempdir(parent)
    try
        cd(fn, tmpdir)
    finally
        rm(tmpdir, recursive=true)
    end
end
