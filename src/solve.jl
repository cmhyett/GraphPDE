include("../types/networkTypes.jl");

# sg is the initial state,
# innerUpdate and boundaryUpdate are both functions that consume a SerialGraph
#  and modify in place du
function createODEProblem(sg_0::SerialGraph, innerUpdate!, boundaryUpdate!, tspan)
    function f!(du, u, p, t)
        innerUpdate!(du, sg_0, p, t);
        boundaryUpdate!(du, sg_0, p, t);
    end

    return ODEProblem(f!, sg_0, tspan);
end
