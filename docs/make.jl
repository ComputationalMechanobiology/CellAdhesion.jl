using Documenter
using Literate
using Revise
using CellAdhesion

DocMeta.setdocmeta!(CellAdhesion, :DocTestSetup, :(using CellAdhesion); recursive = true)

"""
    docprepare()

Convert all source documentation (.jl and .md files)
into markdown, ready for documentation build, and 
relevant .jl files into Jupyter notebooks ready for
deployment to notebooks branch.
"""
function docprepare()
    # remove remnants of previous build and recreate staging dir
    rm("docs/staging-docs", force=true, recursive=true)
    mkdir("docs/staging-docs")

    # copy readme to staging-docs, remove Logo image
    write("docs/staging-docs/index.md", 
            open("README.md") do input
                #readuntil(input, "<!-- delim -->", keep = true)
                read(input)
            end)
    
    # copy assets to staging directory
    cp("docs/src/assets", "docs/staging-docs/assets")

    # iterate through src and convert/copy as appropriate
    for file in readdir("docs/src")
        if endswith(file, "md")
            cp("docs/src/$file", "docs/staging-docs/$file")
        elseif endswith(file, "jl")
            Literate.markdown("docs/src/$file", "docs/staging-docs/"; documenter=true)
        end
    end
end

function notebookprepare()
    # create notebook staging dir
    mkdir("docs/staging-docs/notebooks")

    # copy assets to notebooks staging directory
    cp("docs/src/assets", "docs/staging-docs/notebooks/assets")

    for file in readdir("docs/src")
        if endswith(file, "jl")
            Literate.notebook("docs/src/$file", "docs/staging-docs/notebooks/")
        end
    end
end

function maindocbuilder()
    # prepare doc files from source
    docprepare()

    # prepare notebook files from source
    notebookprepare()

    # build docs from staging area
    makedocs(modules=[CellAdhesion],
            doctest = true, clean=true,
            format = Documenter.HTML(),
            sitename ="CellAdhesion.jl",
            source = "staging-docs",
            authors = "Alessandra Bonfati, Alexandre Kabla",
            pages = ["Home" => "index.md", 
            	     "Example: Single junction" => "example_junction.md",
            	     "Architecture" => "architecture.md",
            	     "API" => "API.md"])

    deploydocs(repo = "github.com/ComputationalMechanobiology/CellAdhesion.jl.git",
               deps = nothing,
               make = nothing,
               target = "build")
               

end

maindocbuilder()
