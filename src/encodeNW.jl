function encodeGraph(g::Graph)
    ndims,_ = size(g.edges[1].initialValues);
    lEdges = sum([length(g.edges[i]) for i in 1:length(g.edges)]);
    lNodes = length(g.nodes);
    u = zeros(n_dims, lEdges + lNodes);
    for 
end
