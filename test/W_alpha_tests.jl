@testset "Proposition 7.10, a_i even" begin
    for a in [[8,2,4,2],[2,2],[8,4,2],[8,4,2],[6,6,10,2],[2,8,8,2],[2,2,2,2,2],[8]]
 #       print('*')
        m = sum(a)
        N = length(a)
        W = QQMatrixBaryNAxis(a)
        # first proof of proposition 7.10
        Q = QQMatrixQ(m)
        CoSqW = transpose(inv(W))*W
        @test CoSqW == -identity_matrix(QQ,m)+QQ(2,N)*Q*ones_matrix(QQ,m,m)
        C = QQ(1,2)*identity_matrix(QQ,m) + QQMatrixU(m)
        CoSqC = transpose(inv(C))*C
        @test CoSqC == -identity_matrix(QQ,m)+2*Q*ones_matrix(QQ,m,m)
        @test jordan_normal_form(CoSqW)[1] == jordan_normal_form(CoSqC)[1]
        # second proof of proposition 7.10
        P1 = identity_matrix(QQ,m) - (transpose(QQMatrixU(m) - block_diagonal_matrix([QQMatrixU(ai) for ai in a])))*QQMatrixQ(m)
        @test P1*W*transpose(P1) == QQ(1,2*N^2)*(identity_matrix(QQ,m) + (N+1)*QQMatrixU(m) - (N-1)*transpose(QQMatrixU(m)))
        P2 = identity_matrix(QQ,m) + (N-1)*sum([QQMatrixE(m,i,j)*QQMatrixQ(m) for j in (1:m) for i in (Int(2*floor(j/2)+1):m)])
        A = QQMatrixSigAxis(m)
        @test P2*P1*W*transpose(P2*P1) == A
    end
end

@testset "Proposition 7.10, a_1 odd" begin
    for a in [[7,2,4,2],[1,2],[3,4,2],[7,3,5],[5,6,10,2],[1,8,7,2],[1,2,2,2,1],
              [7,2,1,2,12],[1,1,1],[3,12,2],[1,1,1,3,5],[3,6,1,2],[3,1,7,2],[1,1,1,1,1]]
#        print('*')
        m = sum(a)
        N = length(a)
        W = QQMatrixBaryNAxis(a)
        # first proof of proposition 7.10
        P1 = identity_matrix(QQ,m) + sum([(-1)^j*QQMatrixE(m,i,j) for i in (a[1]+1:m) for j in (1:a[1])])
        temp = QQ(1,2*N)*block_diagonal_matrix([QQMatrixU(ai) - transpose(QQMatrixU(ai)) for ai in a])
        @test P1*W*transpose(P1) == temp+QQ(1,2*N^2)*block_diagonal_matrix([ones_matrix(QQ,a[1],a[1]),zero_matrix(QQ,m-a[1],m-a[1])])
        Z = QQMatrixU(a[1]) - transpose(QQMatrixU(a[1])) + QQ(1,N)*ones_matrix(QQ,a[1],a[1])
        @test inv(Z) == QQMatrixQ(a[1])*(QQMatrixU(a[1]) - transpose(QQMatrixU(a[1])) + N*ones_matrix(QQ,a[1],a[1]))*QQMatrixQ(a[1])
        C = QQ(1,2)*identity_matrix(QQ,a[1]) + QQMatrixU(a[1])
        CoSqZ = transpose(inv(Z))*Z
        CoSqC = transpose(inv(C))*C
        @test CoSqZ == CoSqC
        # second proof of proposition 7.10
        if a[1] != 1
            P2 = N*identity_matrix(QQ,a[1]) + (N-1)*sum(QQMatrixE(a[1],i,2*j-1)-QQMatrixE(a[1],i,2*j) for j in (1:div(a[1]-1,2)) for i in (2*j:a[1]))
            G = (P1*W*transpose(P1))[1:a[1],1:a[1]]
            @test P2*G*transpose(P2) == C
        end
    end
end





