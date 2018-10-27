push!(LOAD_PATH,"../src/")

using HermiteGF
using Documenter
using Plots # to not capture precompilation output

makedocs(modules=[HermiteGF],
         doctest = false,
         format = :html,
         sitename = "HermiteGF.jl",
         pages = ["Functions"    => "index.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/JuliaVlasov/HermiteGF.jl.git",
    julia  = "1.0",
    osname = "linux"
 )
