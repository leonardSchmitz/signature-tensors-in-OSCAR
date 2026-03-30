using Documenter
using SignatureTensors
using DocumenterCitations

# Bibliografía

bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"), style = :numeric)
# Generar la documentación
makedocs(
    sitename = "SignatureTensors.jl",
    modules = [SignatureTensors],
    pages = [
        "Home" => "index.md",
        "Documentation" => "api.md",
        "References" => "references.md",
    ],
    checkdocs = :warn,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/leonardSchmitz/signature-tensors-in-OSCAR.git",
    branch = "gh-pages",
    devurl = "docs",
    forcepush = true
)

