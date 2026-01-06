# Example 8.2

F, s = free_trunc_sig_alg_multiv(3,2)        # FreeTruncSigAlgMultiv for k=3, N=2

inv(free_sig_from_sample(1,F))               # inverse of s₁⁽³⁾ + s₁⁽²⁾+ s₁⁽¹⁾ + 1    
#-s₁⁽¹⁾^3 + s₁⁽¹⁾^2 + s₁⁽¹⁾*s₁⁽²⁾ + s₁⁽²⁾*s₁⁽¹⁾ - s₁⁽¹⁾ - s₁⁽²⁾ - s₁⁽³⁾ + 1
