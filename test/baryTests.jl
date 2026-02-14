@testset "Barycentric Tests" begin
    # ---------------------------------------------------------
    # Setup inicial
    # ---------------------------------------------------------
    d = 4        # path dimension
    k = 5        # truncation level
    #alpha = [4, 5, 7, 4]
    #N = length(alpha)

    #R, a = polynomial_ring_sig_transform(d, d)
    #TTSd = trunc_tensor_seq(R, k, d)

    #A = [R.(rand(0:100, d, alpha[i])) for i in 1:N]
    #axP = [sig_axis(trunc_tensor_seq(R, k, alpha[i])) for i in 1:N]
    #sX = [matrix_tensorSeq_congruence(A[i], axP[i]) for i in 1:N]

    #sX1, sX2, sX3, sX4 = sX[1:4]
   
    T = TruncatedTensorAlgebra(QQ,d,k);
    ms = rand((1:7),4);
    As = [QQ.(rand((-20:20),d,m)) for m in ms];
    sX = [sig(T,:pwln,coef=A) for A in As];
    sX1, sX2, sX3, sX4 = sX[1:4];

    # ---------------------------------------------------------
    # Theorem 4.3
    # ---------------------------------------------------------
    @testset "Theorem 4.3" begin
        @test sX4 * bary([sX1, sX2]) * sX3 == bary([sX4*sX1*sX3, sX4*sX2*sX3])
        @test bary([sX1, sX2, sX3]) * sX4 == bary([sX1*sX4, sX2*sX4, sX3*sX4])
        @test sX4 * bary([sX1, sX2, sX3]) == bary([sX4*sX1, sX4*sX2, sX4*sX3])
        @test bary([inv(sX1), inv(sX2), inv(sX3), inv(sX4)]) == inv(bary([sX1, sX2, sX3, sX4]))
    end


    # ---------------------------------------------------------
    # Proposition 4.4
    # ---------------------------------------------------------
    @testset "Proposition 4.4" begin
        @test bary([sX1, sX2]) == sX1 * exp(QQ(1,2) * log(inv(sX1) * sX2))
        @test bary([one(T), sX2]) == exp(QQ(1,2) * log(sX2))
        #@test bary_2samples(sX1, sX2) == bary([sX1, sX2])
    end

    # ---------------------------------------------------------
    # Proposition 4.5 / Example 4.6
    # ---------------------------------------------------------
    @testset "Proposition 4.5 / Example 4.6" begin
        @test bary([sX1^2, one(T)]) == sX1
        @test bary([sX1, inv(sX1)]) == one(T)
        @test bary([sX1^6, sX1^2, sX1]) == sX1^3
    end

    # ---------------------------------------------------------
    # Proposition 4.8
    # ---------------------------------------------------------
    @testset "Proposition 4.8" begin
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

    # ---------------------------------------------------------
    # Optional: timing comments
    # ---------------------------------------------------------
    # @time bary_2samples(sX1, sX2)
    # @time bary([sX1, sX2])
end


