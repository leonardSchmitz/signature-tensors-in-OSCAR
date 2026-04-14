export 
    polynomial_path_iis_signature


################################################################################
#  Integration
################################################################################

function exponent_vectors(a::MPolyRingElem{T}; inplace::Bool = false) where T <: RingElement
   return Generic.MPolyExponentVectors(a, inplace=inplace)
end

function exponent_vectors(::Type{Vector{S}}, a::MPolyRingElem{T}; inplace::Bool = false) where {T <: RingElement, S}
   return Generic.MPolyExponentVectors(Vector{S}, a, inplace=inplace)
end

function integration(f::MPolyRingElem{T}, j::Int) where T <: RingElement
   R = parent(f)
   iterz = zip(coefficients(f), exponent_vectors(f))
   Ctx = Generic.MPolyBuildCtx(R)
   for (c, v) in iterz
      if v[j] >= 1
         prod = QQ(1,(v[j]+1))*c
         v[j] += 1
         push_term!(Ctx, prod, v)
      end
   end
   return finish(Ctx)
end

function integration(f::MPolyRingElem{T}, x::MPolyRingElem{T}) where T <: RingElement
   return integration(f, var_index(x))
end


function polynomial_path_iis_signature(A::AbstractMatrix, k::Int; R)

    d, m = size(A)

    T, tvec = polynomial_ring(R, [:t])
    t = tvec[1]

    X = Vector{typeof(t)}(undef, d)
    for i in 1:d
        p = zero(T)
        for j in 1:m
            p += A[i,j] * t^j
        end
        X[i] = p
    end

    DX = [derivative(X[i], t) for i in 1:d]

    results = Vector{Array{typeof(one(R))}}(undef, k+1)
    results[1]= fill(one(R), ())    # array 0-dimensional

    for l in 1:k

        res = Array{typeof(one(R))}(undef, ntuple(_->d, l)...)

        for idx in CartesianIndices(res)

            w = Tuple(idx)

            P = X[w[1]]
            for r in 2:l
                Q = P * DX[w[r]]
                P = integration(Q, t)
            end

            res[idx] = evaluate(P, [one(R)])
        end

        results[l+1] = res
    end

    return results
end




















