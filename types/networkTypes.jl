using RecursiveArrayTools

struct Edge
    length;
    fromNode;
    toNode;
    value; #over spatial discretization size(edge.value)=(ndim,nx)
    params; #properties
end

struct Node
    incidentEdges;
    value; #for single time, size(node.value)=(ndim)
end

mutable struct Graph
    edges::Vector{Edge};
    nodes::Vector{Node};
    graphParams;
    function Graph(edges::Vector{Edge},
                   nodes::Vector{Node},
                   params)
        numEdges = length(edges);
        numNodes = length(nodes);
        ndims_edges = [size(edges[i].value)[1] for i in 1:numEdges];
        ndims_nodes = [size(nodes[i].value)[1] for i in 1:numNodes];

        [@assert ndims_edges[1] == ndims_edges[i] for i in 2:numNodes]
        [@assert ndims_nodes[1] == ndims_nodes[i] for i in 2:numEdges]
        @assert ndims_edges[1] == ndims_nodes[1];
        ndims = ndims_edges[1];

        if (all([typeof(edges[i].value) <: Matrix for i in 1:numEdges]))
            return new(edges, nodes, params);
        elseif (all([typeof(edges[i].value) <: Vector for i in 1:numEdges]))
            newEdges = [];
            for i in 1:numEdges
                arr = zeros(ndims, ceil(Int, edges[i].length/params["dx"]));
                for j in 1:size(arr)[2]
                    arr[:,j] .= edges[i].value[:];
                end
                push!(newEdges, Edge(edges[i].length,
                                     edges[i].fromNode,
                                     edges[i].toNode,
                                     arr,
                                     edges[i].params));
            end
            return new(newEdges, nodes, params);
        else
            @error "In Graph constructor, inhomogenous edge value types, or scalar value type, give me vector or matrix";
        end
    end
end

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
