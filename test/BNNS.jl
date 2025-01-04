const RAND_TYPES = [BFloat16, Float16, Float32, Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64,
                    UInt64]
const RANDN_TYPES = [BFloat16, Float16, Float32]
const INPLACE_TUPLES = [[(rand!, T) for T in RAND_TYPES];
                        [(randn!, T) for T in RANDN_TYPES]]
const OOPLACE_TUPLES = [[(BNNS.rand, rand, T) for T in RAND_TYPES];
                        [(BNNS.randn, rand, T) for T in RANDN_TYPES]]

@testset "random" begin
    # in-place
    @testset "in-place" begin
        rng = BNNS.bnns_rng()

        @testset "$f with $T" for (f, T) in INPLACE_TUPLES
            # d == 2 and d == 3 are to hit the test cases where sizeof(A) <= 4
            @testset "$d" for d in (2, 3, (3, 3), (3, 3, 3), 16, (16, 16), (16, 16, 16), (1000,), (1000,1000))
                A = Array{T}(undef, d)

                # specifie BNNS rng
                fill!(A, T(0))
                f(rng, A)
                @test !iszero(collect(A))
            end

            @testset "0" begin
                A = Array{T}(undef, 0)

                # specified BNNS rng
                fill!(A, T(0))
                f(rng, A)
                @test Array(A) == fill(1, 0)
            end
        end
    end
    # out-of-place
    @testset "out-of-place" begin
        @testset "$fr with implicit type" for (fm, fr, T) in
                                             ((BNNS.rand, Random.rand, Float32), (BNNS.randn, Random.randn, Float32))
            rng = BNNS.bnns_rng()
            @testset "args" for args in ((0,), (1,), (3,), (3, 3), (16,), (16, 16), (1000,), (1000,1000))
                # default_rng
                A = fm(args...)
                @test eltype(A) == T

                # specified MPS rng
                B = fr(rng, args...)
                @test eltype(B) == T
            end

            @testset "scalar" begin
                a = fm()
                @test typeof(a) == T
                b = fr(rng)
                @test typeof(b) == T
            end
        end

        # out-of-place, with type specified
        @testset "$fr with $T" for (fm, fr, T) in OOPLACE_TUPLES
            rng = BNNS.bnns_rng()
            @testset "$args" for args in ((T, 0),
                                          (T, 1),
                                          (T, 3),
                                          (T, 3, 3),
                                          (T, (3, 3)),
                                          (T, 16),
                                          (T, 16, 16),
                                          (T, (16, 16)),
                                          (T, 1000),
                                          (T, 1000, 1000),)
                # default_rng
                A = fm(args...)
                @test eltype(A) == T

                # specified RNG rng
                B = fr(rng, args...)
                @test eltype(B) == T
            end

            @testset "scalar" begin
                a = fm(T)
                @test typeof(a) == T
                b = fr(rng, T)
                @test typeof(b) == T
            end
        end
    end

    ## seeding
    @testset "Seeding" begin
        @testset "$d" for d in (1, 3, (3, 3, 3), 16, (16, 16), (16, 16, 16), (1000,), (1000,1000), (3,3,3,3), (3,3,3,3,3), (3,3,3,3,3,3))
            rng = BNNS.bnns_rng(1)
            a = rand(rng, Float32, d)
            Random.seed!(rng, 1)
            b = rand(rng, Float32, d)
            @test a == b
        end
    end
end # testset
