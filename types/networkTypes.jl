using RecursiveArrayTools

struct Edge
    length;
    fromNode;
    toNode;
    value; #over spatial discretization size(edge.value)=(ndim,nx)
    params; #properties
end

struct Node
    incomingEdges;
    outgoingEdges;
    value; #for single time, size(node.value)=(ndim)
end

#TODO this struct is likely unnecessary..
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

        [@assert ndims_edges[1] == ndims_edges[i] for i in 2:numEdges]
        [@assert ndims_nodes[1] == ndims_nodes[i] for i in 2:numNodes]
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
    numEdges::Int;
    numNodes::Int;
    incomingEdges::Array{Array{Int, 1}, 1}; #length=numNodes
    outgoingEdges::Array{Array{Int, 1}, 1}; #length=numNodes
    edgeProperties::Array{Any, 1}; #length=numEdges
    nodeProperties::Array{Any, 1}; #length=numNodes
end
    
SerialGraph(u::AbstractArray{T,N}, numEdges, numNodes, incomingEdges, outgoingEdges,
            edgeProperties, nodeProperties) where {T,N} =
                SerialGraph{eltype(u),N}(u, numEdges, numNodes, incomingEdges, outgoingEdges,
                                         edgeProperties, nodeProperties);

function SerialGraph(g::Graph)
    return SerialGraph(ArrayPartition([g.edges[i].value for i in 1:length(g.edges)]...,
                                      [g.nodes[i].value for i in 1:length(g.nodes)]...),
                       length(g.edges),
                       length(g.nodes),
                       [g.nodes[i].incomingEdges for i in 1:length(g.nodes)],
                       [g.nodes[i].outgoingEdges for i in 1:length(g.nodes)],
                       [g.edges[i].params for i in 1:length(g.edges)],
                       []);#nodeProperties
end


Base.size(var::SerialGraph) = size(var.u);
Base.getindex(var::SerialGraph, i::Int) = var.u[i];
Base.getindex(var::SerialGraph, I::Vararg{Int,N}) where {N} = var.u[I...];
Base.getindex(var::SerialGraph, ::Colon) = var.u[:];
Base.getindex(var::SerialGraph, kr::AbstractRange) = var.u[kr];

Base.setindex!(var::SerialGraph, v, i::Int) = (var.u[i] = v);
Base.setindex!(var::SerialGraph, v, I::Vararg{Int,N}) where {N} = (var.u[I...] = v);
Base.setindex!(var::SerialGraph, v, ::Colon) = (var.u[:] .= v);
Base.setindex!(var::SerialGraph, v, kr::AbstractRange) = (var.u[kr] .= v);
Base.similar(var::SerialGraph) = SerialGraph(similar(var.u),
                                             var.numEdges,
                                             var.numNodes,
                                             var.incomingEdges,
                                             var.outgoingEdges,
                                             var.edgeProperties,
                                             var.nodeProperties);
Base.similar(var::SerialGraph,::Type{T}) where T = SerialGraph(similar(var.u,T),
                                                               var.numEdges,
                                                               var.numNodes,
                                                               var.incomingEdges,
                                                               var.outgoingEdges,
                                                               var.edgeProperties,
                                                               var.nodeProperties);
