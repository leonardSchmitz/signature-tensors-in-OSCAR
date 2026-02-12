@testset "tensor learning for :iis" begin

    @testset "Pwlin recovery with 1 solution d=m=4 and k=3" begin
        d =4; k = 3;
        R, a = polynomial_ring(QQ, :a => (1:d,1:d));
        C= sig(TruncatedTensorAlgebra(R,d,k),:axis); 
        nr_rec = 3;
        for t in (1:nr_rec)
          A = generic_transform_GL(d);
          G = A*C;
          F = a*C;
          I = ideal(R,vec(F-G));
          @test dim(I) == 0
          @test degree(I) == 1
          @test normal_form.(a, Ref(I)) == Matrix{QQMPolyRingElem}(R.(A))
        end
    end


    @testset "Pwlin recovery with 4 solutions d= 2; k= 4; m = 4" begin
        # where we get the correct via factorization and elimination
        d= 2; k= 4; m = 4; 
        A=QQ[6 -2 6 -10; 7 -4 10 -4];
        R, a = polynomial_ring(QQ, :a => (1:d, 1:m));
        S= sig(TruncatedTensorAlgebra(QQ,d,k),:pwln,coef=A); 
        C= sig(TruncatedTensorAlgebra(QQ,m,k),:axis); 
        I= ideal(R,vec(S-a*C));
        @test dim(I) == 0
        @test degree(I) == 4
        Iel = eliminate(I,vec(a)[1:end-1])[1];
        f = collect(factor(Iel))[1][1];
        @test normal_form.(a, Ref(I+ideal(R,[f]))) == Matrix{QQMPolyRingElem}(R.(A))
    end

end