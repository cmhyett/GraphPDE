using Test;
include("../types/networkTypes.jl");

@testset "Constructor" begin
    # basic network, one edge, two nodes.
    e = Edge(5, 1, 2, zeros(1,5), []);
    n1 = Node([1], zeros(1));
    n2 = Node([1], zeros(1));
    @test typeof(Graph([e], [n1, n2], [])) == Graph

    # wrong number of dimensions (2 != 1)
    e = Edge(5, 1, 2, zeros(2,5), []);
    n = Node([1], zeros(1));
    @test_throws AssertionError Graph([e], [n], [])

    # wrong number of dimensions (1 != 2)
    e = Edge(5, 1, 2, zeros(1,5), []);
    n = Node([1], zeros(2));
    @test_throws AssertionError Graph([e], [n], [])

    # network with specified discretization
    ndims = 3;
    e = Edge(5, 1, 2, zeros(ndims), []);
    n1 = Node([1], zeros(ndims));
    n2 = Node([1], zeros(ndims));
    @test typeof(Graph([e], [n1, n2], Dict("dx"=>0.5))) == Graph
    g = Graph([e], [n1, n2], Dict("dx"=>0.5));
    @test size(g.edges[1].value) == (ndims,5/0.5)
end

