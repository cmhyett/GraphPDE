using Test;

@testset "Constructor" begin
    e = Edge(5, 1, 2, zeros(1,5), []);
    n = Node([1], zeros(1));
    @test typeof(Graph([e], [n], [])) == Graph

    e = Edge(5, 1, 2, zeros(2,5), []);
    n = Node([1], zeros(1));
    @test_throws AssertionError Graph([e], [n], [])

    e = Edge(5, 1, 2, zeros(1,5), []);
    n = Node([1], zeros(2));
    @test_throws AssertionError Graph([e], [n], [])

    e = Edge(5, 1, 2, zeros(5), []);
    n = Node([1], zeros(2));
    @test_throws AssertionError Graph([e], [n], [])

    # lets say I want to initialize a network
    e = Edge(5, 1, 2, zeros(5), []);
    n = Node([1], zeros(5));
    @test_throws ArgumentError Graph([e], [n], []) #no field dx!

    e = Edge(5, 1, 2, zeros(5), []);
    n = Node([1], zeros(5));
    @test typeof(Graph([e], [n], Dict("dx"=>1))) == Graph
    g = Graph([e], [n], Dict("dx"=>1));
    @test size(g.edges[1].value) == (5,5)
end

