module signatureTensors

using Oscar;

greet() = print("Hello World!")

include("typesetting.jl")
include("tensorOperations.jl")
include("FreeTruncSigAlgMultiv.jl")
include("TruncTensorSeq.jl")           # to be removed soon 
include("matrixOperations.jl")
include("TruncatedTensorAlgebra")


end 
