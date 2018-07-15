let

@testset "possiblesums"
    a = [1,1,3,4]
    mina = minimum(a); maxa = maximum(a)
    for n in 1:3
        allsum = all_possible_sums(a, n)
        @test length(allsum) == length(a)^n
        for i in 1:length(a)^n
            @test mina*n ≤ allsum[i][1] ≤ maxa*n
        end
    end
end

end
