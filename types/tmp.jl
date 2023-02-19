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
