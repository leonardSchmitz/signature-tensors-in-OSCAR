@testset "Barycentric Tests" begin
    # CDMSSU24 : https://epubs.siam.org/doi/abs/10.1137/23M159024X
    # AS25 : https://doi.org/10.48550/arXiv.2509.07815
    d = 4        # path dimension
    k = 5        # truncation level
   
    T = TruncatedTensorAlgebra(QQ,d,k);
    ms = rand((1:7),4);
    As = [QQ.(rand((-20:20),d,m)) for m in ms];
    sX = [sig(T,:pwln,coef=A) for A in As];
    sX1, sX2, sX3, sX4 = sX[1:4];

    @testset "Theorem 4.3 in AS25" begin
        @test sX4 * bary([sX1, sX2]) * sX3 == bary([sX4*sX1*sX3, sX4*sX2*sX3])
        @test bary([sX1, sX2, sX3]) * sX4 == bary([sX1*sX4, sX2*sX4, sX3*sX4])
        @test sX4 * bary([sX1, sX2, sX3]) == bary([sX4*sX1, sX4*sX2, sX4*sX3])
        @test bary([inv(sX1), inv(sX2), inv(sX3), inv(sX4)]) == inv(bary([sX1, sX2, sX3, sX4]))
    end

    @testset "Proposition 4.4 in AS25" begin
        @test bary([sX1, sX2]) == sX1 * exp(QQ(1,2) * log(inv(sX1) * sX2))
        @test bary([one(T), sX2]) == exp(QQ(1,2) * log(sX2))
        #@test bary_2samples(sX1, sX2) == bary([sX1, sX2])
    end

    @testset "Proposition 4.5 / Example 4.6 in AS25" begin
        @test bary([sX1^2, one(T)]) == sX1
        @test bary([sX1, inv(sX1)]) == one(T)
        @test bary([sX1^6, sX1^2, sX1]) == sX1^3
    end

    @testset "Proposition 4.8 in AS25" begin
        # i)
        @test bary([sX1, sX1]) == sX1
        @test bary([sX1, sX1, sX1]) == sX1
        @test bary([sX1, sX1, sX1, sX1]) == sX1
        # ii)
        @test bary([sX1, sX2, sX3, sX4]) == bary([sX4, sX1, sX2, sX3])
        @test bary([sX1, sX2, sX3, sX4]) == bary([sX4, sX3, sX2, sX1])
        # iii)
        @test bary_all_but_last_samples_fixed_inverse([sX1, sX2, sX3], bary(sX)) == sX4
    end

    d = 4;       # path dimension
    k = 2;        # truncation level
    T = TruncatedTensorAlgebra(QQ,d,k);
    ms = rand((1:7),4);
    As = [QQ.(rand((-20:20),d,m)) for m in ms];
    sX = [sig(T,:pwln,coef=A) for A in As];
    sX1, sX2, sX3, sX4 = sX[1:4];

    @testset "Remark 4.10 in AS25" begin
        @test bary(sX) == bary(sX,algorithm=:CDMSSU24trunc2)
    end

    @testset "Proposition 4.11 in AS25" begin
        @test bary(sX) == bary(sX,algorithm=:AS25trunc2)
    end

    @testset "Theorem 7.2 in AS25" begin
        for alpha in [[2,2],[3,1,2],[2,2,1,2],[2,1,1,2], 
                      [1,2,1,2],[1,1,1],[1,1,1,3],[1,1,1,1,1]]
            N = length(alpha); sum_alpha = sum(alpha);
            d = sum_alpha;
            odds = length([ai for ai in alpha if isodd(ai)]);
            if odds == sum_alpha
              Bd2 = sum_alpha;
            else 
              Bd2 = sum_alpha - odds + 1;
            end 
            R, a = polynomial_ring(QQ, :a => (1:d,1:Bd2));                          # d = 2, m = Bd2
            TTSd = TruncatedTensorAlgebra(R,d,k);
            A = [R.(rand(0:100, d, alpha[i])) for i in (1:N)];
            sX = [sig(TTSd ,:pwln,coef=A[i]) for i in (1:N)];
            sY = bary(sX);
            s_pwln = sig(TTSd,:pwln,coef=a);          
            I = ideal(R,vec(s_pwln - sY));
            @test is_one(I) != true
        end
    end
end



