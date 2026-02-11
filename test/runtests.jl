using Test
using Oscar
using SignatureTensors   

println("Running all tests in the tests/ folder...")

# Current test folder
test_dir = @__DIR__

# Get all .jl files except runtests.jl
test_files = filter(f -> endswith(f, ".jl") && f != "runtests.jl", readdir(test_dir))

# Run each file inside a global testset
@testset "Running all tests" begin
    for file in sort(test_files)
        println("\n=== Including file: ", file, " ===")
        include(joinpath(test_dir, file))
    end
end

println("\nAll tests completed.")
