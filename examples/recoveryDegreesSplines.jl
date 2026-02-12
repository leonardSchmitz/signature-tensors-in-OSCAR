using SignatureTensors
using Oscar


function ourFilter(v,r)
  canonical(c) = min(c, reverse(c))
  res = unique(canonical, v)
  res = filter(c -> all(x -> x >= r, c), res)
  res = filter(c -> !(any(c[i] == r && c[i+1] == r for i in 1:length(c)-1)), res)
end

function comp_S_Mell(M::Int,ell::Int,r::Int)
  res = [m for m in collect(compositions(M)) if length(m)==ell]
  #return unique(res)
  return ourFilter(res ,r)
end

function ourTable(r::Int)
  # valid M and ell  
  P = lattice_points(polyhedron(([-1 0;0 -1;-2 2*r+1],[-1;-1;2]),([1 -r],[4 - r])))
  #P = lattice_points(polyhedron(([-1 0;0 -1; 0 1],[-1;-1;9-2*r]),([1 -r],[4 - r])))
  #P = filter(c -> c[2]!=0, unique(P)) # no words of length 0 
  #if r >= 1 
  #  P = filter(c -> c[2]!=1, unique(P)) # no words of length 1
  #end 
  res = vcat([comp_S_Mell(Int(x[1]),Int(x[2]),r) for x in P]...)
  return res
end

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
