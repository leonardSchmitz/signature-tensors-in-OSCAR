using Oscar, SignatureTensors

# Define a truncated tensor algebra: dimension d=2, truncation level k=3
d, k = 2, 3
T = TruncatedTensorAlgebra(QQ, d, k)

# Signature of the canonical axis path
C=sig(T, :axis)

# Signature of a polynomial path t ↦ (t + 2t², 3t + 4t²)
S=sig(T, :pwln, coef = QQ[1 2; 3 4])

# Path recovery from a signature tensor
recover(S,Co=C)