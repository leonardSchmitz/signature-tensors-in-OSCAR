
@testset "Matrix-Tensor mode-product" begin
  number_tests = 10;
  ms = rand((1:15),number_tests); # maximal 7 number of segments
  ds = rand((1:15),number_tests); # maximal 7 number of segments
  Gs = [QQ.(rand((-20:20),ms[i],ms[i],ms[i])) for i in (1:number_tests)];
  As = [QQ.(rand((-20:20),ds[i],ms[i])) for i in (1:number_tests)];
  for i in (1:number_tests)
    for mu in (1:3)
      @test matrix_tensor_multiply_3(As[i], Gs[i], mu) == matrix_tensor_multiply(As[i], Gs[i], mu);
      #@test matrix_tensor_multiply_TA(As[i], Gs[i], mu) == matrix_tensor_multiply(As[i], Gs[i], mu);
    end
  end

  As = [generic_transform_GL(ms[i]) for i in (1:number_tests)];
  for i in (1:number_tests)
    for mu in (1:3)
      AG = matrix_tensor_multiply(As[i], Gs[i], mu);
      @test matrix_tensor_multiply(inv(As[i]), AG, mu) == Gs[i];
    end
  end
end