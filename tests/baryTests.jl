@testset "Barycentric Tests" begin
    # ---------------------------------------------------------
    # Setup inicial
    # ---------------------------------------------------------
    d = 4        # path dimension
    k = 5        # truncation level
    alpha = [4, 5, 7, 4]
    N = length(alpha)

    R, a = polynomial_ring_sig_transform(d, d)
    TTSd = trunc_tensor_seq(R, k, d)

    A = [R.(rand(0:100, d, alpha[i])) for i in 1:N]
    axP = [sig_axis(trunc_tensor_seq(R, k, alpha[i])) for i in 1:N]
    sX = [matrix_tensorSeq_congruence(A[i], axP[i]) for i in 1:N]

    sX1, sX2, sX3, sX4 = sX[1:4]

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
    # Basic arithmetic tests
    # ---------------------------------------------------------
    @testset "Basic arithmetic" begin
        @test exp(log(sX1)) == sX1
        @test inv(inv(sX1)) == sX1
        @test sX1 * inv(sX1) == one(TTSd)
        @test inv(sX1) * sX1 == one(TTSd)
        @test sX1^3 == sX1*sX1*sX1
        @test exp(log(sX1) + log(sX1)) == sX1^2
    end

    # ---------------------------------------------------------
    # Proposition 4.4
    # ---------------------------------------------------------
    @testset "Proposition 4.4" begin
        @test bary([sX1, sX2]) == sX1 * exp(QQ(1,2) * log(inv(sX1) * sX2))
        @test bary([one(TTSd), sX2]) == exp(QQ(1,2) * log(sX2))
        @test bary_2samples(sX1, sX2) == bary([sX1, sX2])
    end

    # ---------------------------------------------------------
    # Proposition 4.5 / Example 4.6
    # ---------------------------------------------------------
    @testset "Proposition 4.5 / Example 4.6" begin
        @test bary([sX1^2, one(TTSd)]) == sX1
        @test bary([sX1, inv(sX1)]) == one(TTSd)
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
        @test bary_Nminus1_samples_fixed_inverse([sX1, sX2, sX3], bary(sX)) == sX4
    end

    # ---------------------------------------------------------
    # Optional: timing comments
    # ---------------------------------------------------------
    # @time bary_2samples(sX1, sX2)
    # @time bary([sX1, sX2])
end


