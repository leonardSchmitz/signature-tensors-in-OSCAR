# Example 3.7

F, s = free_trunc_sig_alg_multiv(3,2)        # FreeTruncSigAlgMultiv for k=3, N=2

z1 = free_sig_from_sample(1,F)
z2 = free_sig_from_sample(2,F)

z1*z2
# s₁⁽¹⁾*s₂⁽¹⁾ + s₁⁽¹⁾*s₂⁽²⁾ + s₁⁽²⁾*s₂⁽¹⁾ + s₁⁽¹⁾ + s₁⁽²⁾ + s₁⁽³⁾ + s₂⁽¹⁾ + s₂⁽²⁾ + s₂⁽³⁾ + 1
