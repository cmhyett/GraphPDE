using Test, DifferentialEquations;
include("../types/networkTypes.jl");

@testset "Constructor" begin
    # basic network, one edge, two nodes.
    e = Edge(5, 1, 2, zeros(1,5), []);
    n1 = Node([1], [], zeros(1));
    n2 = Node([1], [], zeros(1));
    @test typeof(Graph([e], [n1, n2], [])) == Graph

    # wrong number of dimensions (2 != 1)
    e = Edge(5, 1, 2, zeros(2,5), []);
    n = Node([1], [], zeros(1));
    @test_throws AssertionError Graph([e], [n], [])

    # wrong number of dimensions (1 != 2)
    e = Edge(5, 1, 2, zeros(1,5), []);
    n = Node([1], [], zeros(2));
    @test_throws AssertionError Graph([e], [n], [])

    # test network 1
    ndims = 3;
    tspan = (0.0, 1.0)
    n1 = Node([], [1], zeros(ndims));
    e1 = Edge(5, 1, 2, zeros(ndims), []);
    n2 = Node([1], [2,4], zeros(ndims));
    g = Graph([e1], [n1, n2], Dict("dx"=>0.5));
    sg = SerialGraph(g);
    function rhs_test2(df, f, p, t)
        df.u .= f.u;
    end
    prob  = ODEProblem(rhs_test2, sg, tspan)
    sol   = solve(prob, RK4())
    tests = vcat((sol .== 0.0)...);
    for i in 1:length(tests)
        @test tests[i];
    end
end

