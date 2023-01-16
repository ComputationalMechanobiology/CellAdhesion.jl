using Documenter
using CellAdhesion

makedocs(
    sitename = "CellAdhesion",
    format = Documenter.HTML(),
    modules = [CellAdhesion]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/alebonfanti/CellAdhesion.jl.git"
)
