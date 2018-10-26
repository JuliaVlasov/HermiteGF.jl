push!(LOAD_PATH,"../src/")

using HermiteGF
using Documenter
using Plots # to not capture precompilation output

makedocs(modules=[HermiteGF],
         doctest = false,
         format = :html,
         sitename = "HermiteGF.jl",
         pages = ["Introduction"    => "index.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/pnavaro/Splittings.jl.git",
    julia  = "1.0",
    osname = "osx"
 )
