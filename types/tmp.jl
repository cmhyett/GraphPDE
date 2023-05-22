using DifferentialEquations

struct SerialGraph{A,N} <: AbstractArray{A,N}
    u :: ArrayPartition{A}
    edgeInds;
    nodeInds;
end
    
SerialGraph(u::AbstractArray{T,N}, edgeInds, nodeInds) where {T,N} =
    SerialGraph{eltype(u),N}(u, edgeInds, nodeInds);

Base.size(var::SerialGraph) = size(var.u);

Base.getindex(var::SerialGraph, i::Int) = var.u[i];
Base.getindex(var::SerialGraph, I::Vararg{Int,N}) where {N} = var.u[I...];
Base.getindex(var::SerialGraph, ::Colon) = var.u[:];
Base.getindex(var::SerialGraph, kr::AbstractRange) = var.u[kr];

Base.setindex!(var::SerialGraph, v, i::Int) = (var.u[i] = v);
Base.setindex!(var::SerialGraph, v, I::Vararg{Int,N}) where {N} = (var.u[I...] = v);
Base.setindex!(var::SerialGraph, v, ::Colon) = (var.u[:] .= v);
Base.setindex!(var::SerialGraph, v, kr::AbstractRange) = (var.u[kr] .= v);

Base.similar(var::SerialGraph) = SerialGraph(similar(var.u), var.edgeInds, var.nodeInds);
Base.similar(var::SerialGraph,::Type{T}) where T = SerialGraph(similar(var.u,T),var.edgeInds, var.nodeInds);

function rhs_test(f, p, t)
    f
end

xmin   = -2.0*pi
xmax   =  2.0*pi
xnodes =  36000
hx     = (xmax - xmin) / xnodes

xx = range(xmin, stop=xmax, length=xnodes)

x0   = 0
w    = 0.4
A    = 1

f0    = A * exp.( -((xx .- x0) ./ w).^2 )
foo   = SerialGraph(ArrayPartition(f0, f0/10), [1,2], [2,3])

tspan = (0.0, 1.0)


function rhs_test2(df, f, p, t)
    df.u .= f.u;
end
prob  = ODEProblem(rhs_test2, foo, tspan)
sol   = solve(prob, RK4())


                 
# test network
#           3
#        e2/|
#    e1   / |
# 1 ---- 2  |e3
#         \ |
#        e4\|
#           4
ndims = 3;
n1 = Node([1], zeros(ndims));
e1 = Edge(5, 1, 2, zeros(ndims), []);
n2 = Node([1,2,4], zeros(ndims));
e2 = Edge(3, 2, 3, zeros(ndims), []);
n3 = Node([2,3], zeros(ndims));
e3 = Edge(3, 4, 3, zeros(ndims), []);
n4 = Node([3,4], zeros(ndims));
e4 = Edge(3, 4, 2, zeros(ndims), []);
g = Graph([e1, e2, e3, e4], [n1, n2, n3, n4], Dict("dx"=>0.5));
