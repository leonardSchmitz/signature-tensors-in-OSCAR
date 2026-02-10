@testset "Truncated Tensor Algebra Tests" begin
    d = 6        # path dimension
    k = 5        # truncation level
    T = TruncatedTensorAlgebra(QQ,d,k)
    m = 7        # number of segments
    A = QQ.(rand((-20:20),d,m))
    S = sig(T,:pwln,coef=A); 
    C = sig(TruncatedTensorAlgebra(QQ,m,k),:axis); 

    @testset "Constructor TTA" begin
       @test T == TruncatedTensorAlgebra(QQ,d,k,sequence_type=:iis)
       @test sequence_type(T) == :iis
       @test base_dimension(T) == d
       @test base_ring(T) == QQ
       @test truncation_level(T) == k
       @test parent(S) == T
    end

    @testset "Constructor of elements in TTA" begin
       @test zero(T) + zero(T) == zero(T)
       @test zero(T) + S == zero(T)
       @test one(T)*S == S
       @test S*one(T) == S
       @test A*C == S
       @test inv(A)*S == C
       @test sig(T,:axis,algorithm=:AFS19) == sig(T,:axis,algorithm=:Chen) 
       @test sig(T,:pwln,coef=A,algorithm=:congruence) == sig(T,:axis,coef=A,algorithm=:Chen) 
    end
    
    @testset "Chen" begin
        for m1 in 1:m-1
            A1 = A[:,1:m1]
            A2 = A[:,m1+1:end]
            @test hcat(A1,A2) == A
            @test sig(T,:pwln,coef=A) == sig(T,:pwln,coef=A1)*sig(T,:pwln,coef=A2)
        end
    end
end