# Example 8.1

k = 2; d = 2; m = 2
R, a = polynomial_ring_sig_transform(d,m)     # polynomial ring QQ[a₁₁,a₂₁,a₁₂,a₂₂]
TTSm = trunc_tensor_seq(R,k,m)                # creates a TruncTensorSeq

matrix_tensorSeq_congruence(a,sig_axis(TTSm))       # a TruncTensorSeqElem
# 0-dimensional Array{QQMPolyRingElem, 0}:
# 1
# ⊕
# 2-element Vector{QQMPolyRingElem}:
#  a₁₁ + a₁₂
#  a₂₁ + a₂₂
# ⊕
# 2×2 Matrix{QQMPolyRingElem}:
#  1//2*a₁₁^2 + a₁₁*a₁₂ + 1//2*a₁₂^2      1//2*a₁₁*a₂₁ + a₁₁*a₂₂ + 1//2*a₁₂*a₂₂
#  1//2*a₁₁*a₂₁ + a₂₁*a₁₂ + 1//2*a₁₂*a₂₂  1//2*a₂₁^2 + a₂₁*a₂₂ + 1//2*a₂₂^2

