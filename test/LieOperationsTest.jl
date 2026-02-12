@testset "Lie Operation Tests" begin
   
   for k in (2:6)
     for d in (2:6)
       @test dim_free_nil_lie_alg(d,k) == length(generate_lyndon_words(k,Vector((1:d))))
     end 
   end 

end