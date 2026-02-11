
@testset "Truncated Tensor Algebra Tests :iis" begin
    d = 6        # path dimension
    k = 5        # truncation level
    T = TruncatedTensorAlgebra(QQ,d,k)

    @testset "Constructor TTA" begin
       @test T == TruncatedTensorAlgebra(QQ,d,k,sequence_type=:iis)
       @test sequence_type(T) == :iis
       @test base_dimension(T) == d
       @test base_algebra(T) == QQ
       @test truncation_level(T) == k
    end

    function axis_core_3tensor_QQ(_d)
        C = zeros(QQ,_d,_d,_d);
        for al in (1:_d)
          for be in (1:_d)
            for ga in (1:_d)
              if al == be && be == ga
                C[al,be,ga] = QQ(1,6)
              end
              if (al < be && be == ga)||(al == be && be < ga)
                C[al,be,ga] = QQ(1,2)
              end
              if al < be && be < ga
                C[al,be,ga] = one(QQ)
              end
            end
          end
        end
        return C
    end

    @testset "Axis constructor in TTA for QQ" begin
        Caxis_d = sig(T,:axis);
        @test parent(Caxis_d) == T
        @test zero(T) + zero(T) == zero(T)
        @test Caxis_d == sig(T,:axis,algorithm=:Chen)
        @test Caxis_d == sig(T,:axis,algorithm=:AFS19) 
        for i in (2:d-1)
            @test Caxis_d[i] == one(QQ)
            @test Caxis_d[i,i-1] == zero(QQ)
            @test Caxis_d[i,i] == QQ(1,2) 
            @test Caxis_d[i,i+1] == one(QQ)
        end 
        @test axis_core_3tensor_QQ(d) == Caxis_d[:,:,:]
        @test zero(T) + Caxis_d == Caxis_d
        @test one(T)*Caxis_d == Caxis_d
        @test Caxis_d*one(T) == Caxis_d
        @test inv(Caxis_d)*Caxis_d == one(T)
        @test Caxis_d*inv(Caxis_d) == one(T)
        @test inv(inv(Caxis_d)) == Caxis_d
    end

    @testset "Pwln constructor in TTA" begin
        ms = rand((1:7),5);
        As = [QQ.(rand((-20:20),d,m)) for m in ms];
        for A in As
          S = sig(T,:pwln,coef=A); 
          @test S == sig(T,:pwln,coef=A,algorithm=:Chen); 
          @test S == sig(T,:pwln,coef=A,algorithm=:congruence); 
          @test zero(T) + S == S
          @test one(T)*S == S
          @test S*one(T) == S
          @test inv(S)*S == one(T)
          @test S*inv(S) == one(T)
          @test inv(inv(S)) == S
          m = size(A)[2]
          for m1 in 1:m-1     # decompose via Chen
            A1 = A[:,1:m1]
            A2 = A[:,m1+1:end]
            @test hcat(A1,A2) == A
            @test sig(T,:pwln,coef=A) == sig(T,:pwln,coef=A1)*sig(T,:pwln,coef=A2)
          end
        end 
    end        
              
    

    #m = 7        # number of segments
    #A = QQ.(rand((-20:20),d,m)) 
    #C = sig(TruncatedTensorAlgebra(QQ,m,k),:axis); 

    #@test A*Caxis_d == S
    #    Cmono_d = sig(T,:mono)
    #    for Si in [Caxis_d,S,one(T)]
    #        @test zero(T) + Si == zero(T)
    #        @test one(T)*Si == Si
    #        @test Si*one(T) == Si
    #        @test inv(Si)*Si == one(T)
    #        @test Si*inv(Si) == one(T)
    #        @test inv(inv(Si)) == Si
    #    end
    #    for m1 in 1:m-1     # decompose via Chen
    #        A1 = A[:,1:m1]
    #        A2 = A[:,m1+1:end]
    #        @test hcat(A1,A2) == A
    #        @test sig(T,:pwln,coef=A) == sig(T,:pwln,coef=A1)*sig(T,:pwln,coef=A2)
    #    end
    
end