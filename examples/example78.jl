a = [4,6,2]
N = length(a)

B1 = QQMatrixU(a[1]) - transpose(QQMatrixU(a[1]))
B2 = QQMatrixU(a[2]) - transpose(QQMatrixU(a[2]))
B3 = QQMatrixU(a[3]) - transpose(QQMatrixU(a[3]))
B = QQ(1,2*N)*block_diagonal_matrix([B1,B2,B3]) + QQ(1,2*N*N)*ones_matrix(QQ,sum(a),sum(a))

@assert QQMatrixBaryNAxis(a) == B


# show(stdout, "text/latex",block_diagonal_matrix([B1,B2,B3]))
# show(stdout, "text/latex",ones_matrix(QQ,sum(a),sum(a)))
# show(stdout, "text/latex",2*N^2*QQMatrixBaryNAxis(a))

