using Documenter
using SignatureTensors 

makedocs(
    sitename = "SignatureTensors.jl",
    modules = [SignatureTensors],
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
    strict = false,  
)


deploydocs(
    repo = "github.com/leonardSchmitz/signature-tensors-in-OSCAR.git",
)

