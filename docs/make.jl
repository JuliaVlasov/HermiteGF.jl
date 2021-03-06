push!(LOAD_PATH,"../src/")

using HermiteGF
using Documenter

makedocs(modules=[HermiteGF],
         doctest = false,
         format = :html,
         sitename = "HermiteGF.jl",
         pages = ["Documentation" => "index.md",
                  "Functions" => "functions.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/JuliaVlasov/HermiteGF.jl.git",
 )
