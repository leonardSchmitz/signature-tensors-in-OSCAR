@testset "Barycentric Tests" begin
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

end


