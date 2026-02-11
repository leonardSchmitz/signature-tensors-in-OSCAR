@testset "tensor learning for :iis" begin

    @testset "Pwlin recovery with 4 solutions d= 2; k= 4; m = 4" begin
        # where we get the correct via factorization and elimination
        d= 2; k= 4; m = 4; 
        A=QQ[6 -2 6 -10; 7 -4 10 -4];
        R, a = polynomial_ring(QQ, :a => (1:d, 1:m));
        S= sig(TruncatedTensorAlgebra(QQ,d,k),:pwlin,coef=A); 
        C= sig(TruncatedTensorAlgebra(QQ,m,k),:axis); 
        I= ideal(R,vec(S-a*C));
        @test dim(I) == 0
        @test degree == 4
        Iel = eliminate(I,vec(a)[1:end-1])[1];
        f = collect(Iel)[1][1];
        @test normal_form.(a, Ref(I+ideal(R,[f]))) == Matrix{QQMPolyRingElem}(R.(A))
    end

end