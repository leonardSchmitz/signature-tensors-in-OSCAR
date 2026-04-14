"""
Welcome to SignatureTensors v1.0.0

This package allows computing and working with signatures of rough paths
using algebraic methods within OSCAR.

For more information, visit:
https://leonardschmitz.github.io/SignatureTensors.jl/
"""
module SignatureTensors
using Oscar;



#greet() = print("Hello World!")

include("typesetting.jl")
include("LieOperations.jl")
include("tensorOperations.jl")
include("FreeTruncSigAlgMultiv.jl")
include("matrixOperations.jl")
include("polynomialTransformations.jl")
include("TruncatedTensorAlgebra.jl")
include("tensorLearningSch25.jl")


end 
