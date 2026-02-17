using SignatureTensors
using Oscar

# See ALS26 : https://doi.org/10.48550/arXiv.2602.13011

function ourFilter(v,r)
  canonical(c) = min(c, reverse(c))
  res = unique(canonical, v)
  res = filter(c -> all(x -> x >= r, c), res)
  res = filter(c -> !(any(c[i] == r && c[i+1] == r for i in 1:length(c)-1)), res)
end

function comp_S_Mell(M::Int,ell::Int,r::Int)
  res = [m for m in collect(compositions(M)) if length(m)==ell]
  return ourFilter(res ,r)
end

function ourTable(r::Int)
  P = lattice_points(polyhedron(([-1 0;0 -1;-2 2*r+1],[-1;-1;2]),([1 -r],[4 - r])))
  res = vcat([comp_S_Mell(Int(x[1]),Int(x[2]),r) for x in P]...)
  return res
end


# Table 3 and 4 in ALS26 
d = 2; 
k = 4; 
for r in (0:4)
  for m in ourTable(r)
      vm = Vector(m)
      l = length(vm);
      n = sum(m) - (l-1)*r
      R, a = polynomial_ring(QQ, :a => (1:d, 1:n));
      T = TruncatedTensorAlgebra(R,n,k)
      C = sig(T,:pwmon,composition=vm,regularity=r);
      A = generic_transform(d,n);
      I = ideal(R,vec(A*C-a*C));
      LI = leading_monomial.(groebner_basis_f4(I));
      dimI = dim(ideal(R,LI));
      degI = degree(ideal(R,LI));
      println("r=",r,", m=", vm, ", dim=", dimI, ", deg=",degI)
  end
end


function polynomial_ring_sig_transform_geom_splines(_dim::Int,m::Vector{Int},r::Int)
   n = sum(m) - (length(m)-1)*r
   return  polynomial_ring(QQ, :a => (1:_dim,1:n), :rho => (1:(length(m)-1),1:r))
end

function nextBlock_geom_splines(rho,A, r::Int, mi)
    nrows, ncols = size(A)
    R = parent(rho[1,1])
    vecs = zeros(R, ncols, r)
    for i in 1:r
        for j in 1:ncols
            vecs[j, i] = binomial(j, i)*rho[mi,i]
        end
    end
    return matrix(R,A) * matrix(R,vecs)
end

function coreSplineTrafo_geom_splines(rho,m::Vector{Int}, r::Int)
    R = parent(rho[1,1])
    total_dim = sum(m)
    if r == 0
        return identity_matrix(R,total_dim)
    end
    zz = 1
    crb = identity_matrix(R,m[1])
    B = crb
    while zz < length(m)
        nr = size(B, 2)
        start_col = nr - m[zz] + 1
        end_col   = nr
        subB = B[:, start_col:end_col]
        nb = nextBlock_geom_splines(rho,subB, r, zz)
        idpart = identity_matrix(R,m[zz+1] - r)
        B = block_diagonal_matrix([hcat(B, nb), idpart])
        zz += 1
    end
    return B
end


function get_recovery_degree_geometric_spline(d::Int, k::Int, r::Int, m)
  R, a, rho = polynomial_ring_sig_transform_geom_splines(d,m,r)
  T = TruncatedTensorAlgebra(R,sum(m),k)
  C = sig(T,:pwmon,composition=m);
  B = Array(coreSplineTrafo_geom_splines(rho,m,r));
  n = size(B)[1]
  A = generic_transform(d,n);
  v =  vec(vcat(vec(a),vec(QQ.(rand(-10:10, length(m)-1, r)))))
  eval_rho = f -> evaluate(f,v)
  B_eval = eval_rho.(B)
  I = ideal(R,vec(A*(B_eval*C)-a*(B*C)));
  LI = leading_monomial.(groebner_basis_f4(I));
  @assert dim(ideal(R,LI)) == 0;
  return degree(ideal(R,LI));
end

# Table 5 in ALS26
get_recovery_degree_geometric_spline(2,4,1,[2,2,1]) # 32
get_recovery_degree_geometric_spline(2,4,1,[2,1,2]) # 32
get_recovery_degree_geometric_spline(2,4,1,[1,3,1]) # 96
get_recovery_degree_geometric_spline(2,4,2,[3,2]) # 116

# Table 6 in ALS26
get_recovery_degree_geometric_spline(3,3,1,[3,2,1]) # 144
get_recovery_degree_geometric_spline(3,3,1,[3,1,2]) # 84 
get_recovery_degree_geometric_spline(3,3,1,[2,2,2]) # 90
get_recovery_degree_geometric_spline(3,3,2,[4,2]) # 312
get_recovery_degree_geometric_spline(3,3,2,[3,3]) # 168





function get_dim_and_degree_parametric_spline_signature_image_k4(d::Int,r::Int,m)
  l = length(m); k = 4;
  n = sum(m) - r*(l-1);    # dim for core tensor of spline
  R,a,s4,s3,s2,s1 = polynomial_ring(QQ, :a => (1:d, 1:n),:s4 => (1:d,1:d,1:d,1:d),:s3 => (1:d,1:d,1:d),:s2 => (1:d,1:d),:s1 => (1:d));
  #R,a,s1,s2,s3,s4 = polynomial_ring(QQ, :a => (1:d, 1:n),:s1 => (1:d),:s2 => (1:d,1:d),:s3 => (1:d,1:d,1:d),:s4 => (1:d,1:d,1:d,1:d));
  TTSn = TruncatedTensorAlgebra(R,n,k);
  C = sig(TTSn,:pwmon,composition=m,regularity=r);
  aC = a*C;
  s0 = tensor_sequence(one(TTSn))[1];
  s = TruncatedTensorAlgebraElem(parent(aC),[s0,s1,s2,s3,s4]);
  I = ideal(R,vec(aC-s));
  LI = leading_monomial.(groebner_basis_f4(I,eliminate=n*d));
  dimI = dim(ideal(R,LI))-n*d
  degI = degree(ideal(R,LI));
  return dimI,degI
end

# dimension of the image 
# table 1 in ALS26
get_dim_and_degree_parametric_spline_signature_image_k4(2,0,[1,1,1]) # (6, 276) 
get_dim_and_degree_parametric_spline_signature_image_k4(2,1,[1,1,1]) # (2, 16)
get_dim_and_degree_parametric_spline_signature_image_k4(2,1,[2,1]) # (4, 96)
get_dim_and_degree_parametric_spline_signature_image_k4(2,1,[2,1,1]) # (4, 96)
get_dim_and_degree_parametric_spline_signature_image_k4(2,2,[2,2]) # (4, 96)