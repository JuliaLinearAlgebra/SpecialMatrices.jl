using SpecialMatrices
using Documenter
using Literate

# based on:
# https://github.com/jw3126/UnitfulRecipes.jl/blob/master/docs/make.jl

# https://juliadocs.github.io/Documenter.jl/stable/man/syntax/#@example-block
ENV["GKSwstype"] = "100"
ENV["GKS_ENCODING"] = "utf-8"

# generate examples using Literate
lit = joinpath(@__DIR__, "lit")
src = joinpath(@__DIR__, "src")
notebooks = joinpath(src, "notebooks")


DocMeta.setdocmeta!(SpecialMatrices, :DocTestSetup, :(using SpecialMatrices); recursive=true)

execute = true # Set to true for executing notebooks and documenter!
nb = false # Set to true to generate the notebooks
for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath, documenter = execute)
    nb && Literate.notebook(ipath, notebooks, execute = execute)
end


# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages(folder) =
    [joinpath(folder, f) for f in readdir(joinpath(src, folder)) if ismd(f)]

isci = get(ENV, "CI", nothing) == "true"

format = Documenter.HTML(;
    prettyurls = isci,
    edit_link = "master",
    canonical = "https://JuliaMatrices.github.io/SpecialMatrices.jl/stable/",
#   assets = String[],
)

makedocs(;
    modules = [SpecialMatrices],
    sitename = "SpecialMatrices.jl",
    format,
    pages = [
        "Home" => "index.md",
        "Methods" => "methods.md",
        "Examples" => pages("examples")
    ],
)

if isci
    deploydocs(;
        repo = "github.com/JuliaMatrices/SpecialMatrices.jl",
        devbranch = "master",
        devurl = "dev",
        versions = ["stable" => "v^", "dev" => "dev"],
        forcepush = true,
#       push_preview = true,
        # see https://JuliaMatrices.github.io/SpecialMatrices.jl/previews/PR##
    )
end
