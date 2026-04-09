using SignatureTensors
using BenchmarkTools
using Oscar
import Oscar: Oscar, gens, parent, coeff, exponent_vector, monomials, MPolyIdeal
import HomotopyContinuation: HomotopyContinuation, Variable
import LinearAlgebra: LinearAlgebra


# For the following 3 functions, see also https://github.com/taboege/OscarHomotopyContinuation/tree/main
function poly_to_expr(f)
    v = map(x -> Variable(string(x)), gens(parent(f)))
    poly = map(m -> [coeff(f, m), exponent_vector(m, 1)], monomials(f))
    +([*([Rational(c), [v[i]^e for (i,e) in enumerate(a)]...]...) for (c,a) in poly ]...)
end

function System(I::MPolyIdeal; args...)
    HomotopyContinuation.System(map(f -> poly_to_expr(f), gens(I)); args...)
end

function nsolve(I::MPolyIdeal; show_progress=false, args...)
    HomotopyContinuation.solve(System(I); show_progress, args...)
end


d = 2
T = TruncatedTensorAlgebra(QQ,d,3)
A = generic_transform_GL(d);
G = sig(T,:pwln,coef=A);
R, a = polynomial_ring(QQ, :a => (1:d, 1:d));
S = sig(T,:pwln,coef=A);
C = sig(T,:axis);

@benchmark nsolve(ideal(R,vec(S-a*C)[2:end]))
#BenchmarkTools.Trial: 81 samples with 1 evaluation per sample.
# Range (min … max):  52.932 ms … 85.397 ms  ┊ GC (min … max): 0.00% … 0.00%
# Time  (median):     59.044 ms              ┊ GC (median):    0.00%
# Time  (mean ± σ):   61.763 ms ±  7.490 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%
#
#    ▁ ▆▁▃▆▃█ ▃  ▁▃    ▁ ▃  ▁                                   
#  ▇▄█▇██████▇█▄▄██▄▄▄▄█▁█▁▁█▁▇▇▄▇▁▁▄▁▄▁▁▁▁▁▄▄▁▁▁▁▁▁▄▄▁▄▁▁▁▁▁▇ ▁
#  52.9 ms         Histogram: frequency by time        83.1 ms <
#
# Memory estimate: 4.40 MiB, allocs estimate: 104588.

@benchmark groebner_basis(ideal(R,vec(S-a*C)[2:end]))
#BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
# Range (min … max):   81.250 μs … 321.003 ms  ┊ GC (min … max):  0.00% … 39.36%
# Time  (median):      90.375 μs               ┊ GC (median):     0.00%
# Time  (mean ± σ):   122.924 μs ±   3.209 ms  ┊ GC (mean ± σ):  10.28% ±  0.39%
#
#               ▂▅▄▆▅█▇█▅▄                                        
#  ▂▂▂▃▃▃▃▃▃▃▄▅█████████████▆▆▅▅▄▃▃▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂ ▄
#  81.2 μs          Histogram: frequency by time          110 μs <
#
# Memory estimate: 94.70 KiB, allocs estimate: 3386.


d = 3
T = TruncatedTensorAlgebra(QQ,d,3)
A = generic_transform_GL(d);
G = sig(T,:pwln,coef=A);
R, a = polynomial_ring(QQ, :a => (1:d, 1:d));
S = sig(T,:pwln,coef=A);
C = sig(T,:axis);

@time nsolve(ideal(R,vec(S-a*C)[2:end]));
#  743.648308 seconds (100.13 M allocations: 4.112 GiB, 0.17% gc time, 3.08% compilation time)
@time groebner_basis(ideal(R,vec(S-a*C)[2:end]));
#  0.001631 seconds (18.86 k allocations: 500.484 KiB)
