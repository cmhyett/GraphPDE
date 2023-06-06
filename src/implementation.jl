using DifferentialEquations, StatsBase, LinearAlgebra, QuadGK, Plots, Printf;

#=
Single pipe example to show sensitivity capability...
1 is supply, 2 is withdrawal,
1 --> 2
=#
function run()
    L = 10 #km
    dx = 0.05 #km
    x = 0.0:dx:L
    dt = dx;
    ndims = 2;
    tspan = (0.0, 5.0);
    tstops = tspan[1]:dt:tspan[end]
    rho_0 = 100.0;
    phi_0 = 0.0;
    n1 = Node([], [1], [rho_0; phi_0]);
    e1 = Edge(L, 1, 2, zeros(ndims), []);
    n2 = Node([1], [], [rho_0; phi_0]);
    g = Graph([e1], [n1, n2], Dict("dx"=>dx));
    sg = SerialGraph(g);

    sg.u.x[1][1,:] .= rho_0;
    sg.u.x[1][2,:] .= phi_0;
    a = dx/dt;

    A = [0.0 1.0
         a^2 0.0];
    function F(u)
        return A*u;
    end

    function S(u)
        return 0;
    end

    nu = a;
    function LaxFriedrichsFlux(u,i)
        return 0.5 * (F(u[:,i]) + F(u[:,i+1])) - (nu/2)*(u[:,i+1]-u[:,i]);    
    end

    function rightFluxApprox(u,i)
        return LaxFriedrichsFlux(u,i);
    end
    function leftFluxApprox(u,i)
        return LaxFriedrichsFlux(u,i-1);
    end

    function handleBCs!(du,u,p,t)
        # handle edges of interior
        continuation = [u.u.x[2] u.u.x[1][:,1] u.u.x[1][:,2]];
        du.u.x[1][:,1] .= -(1/dx)*(rightFluxApprox(continuation,2)-leftFluxApprox(continuation,2));

        continuation = [u.u.x[1][:,end-1] u.u.x[1][:,end] u.u.x[3]];
        du.u.x[1][:,end] .= -(1/dx)*(rightFluxApprox(continuation,2)-leftFluxApprox(continuation,2));

        # then handle boundary nodes
        continuation = [[u.u.x[2][1]; p[1]] u.u.x[2] u.u.x[1][:,1]];
        du.u.x[2] .= -(1/dx)*(rightFluxApprox(continuation,2)-leftFluxApprox(continuation,2));

        continuation = [u.u.x[1][:,end] u.u.x[3] [u.u.x[3][1]; p[2]]];
        du.u.x[3] .= -(1/dx)*(rightFluxApprox(continuation,2)-leftFluxApprox(continuation,2));
    end

    function dudt!(du,u,p,t)
        for i in 1:u.numEdges
            edge_update = du.u.x[i];
            edge_state = u.u.x[i];
            for j in 2:(size(edge_state)[end])-1
                edge_update[:,j] .= -(1/dx)*(rightFluxApprox(edge_state,j)-leftFluxApprox(edge_state,j));
            end
        end
        handleBCs!(du,u,p,t);
    end

    #params = [0.05; 0.05];
    prob = ODEProblem(dudt!, sg, tspan);
#    sol = Array(solve(prob, RK4(), tstops = tstops, adaptive=false, sensealg=InterpolatingAdjoint()));
    
    # function makePlot(tstep)
    #     densityPlot = plot();
    #     fluxPlot = plot();
    #     label = @sprintf "t = %2.2f" sol.t[tstep]
    #     plot!(densityPlot, midpoints(x), sol.u[tstep].u.x[1][1,:], label=label,
    #           title="Density", ylims=(99.9, 100.1));
    #     plot!(fluxPlot, midpoints(x), sol.u[tstep].u.x[1][2,:], label=label,
    #           title="Flux", ylims=(0.0, 0.1));
    #     return plot(densityPlot, fluxPlot)
    # end

    # anim = @animate for tstep in 1:length(sol.t)
    #     makePlot(tstep)
    # end
    # gif(anim, "anim_fps5.gif", fps=10);
    return prob, tstops
end

prob, tstops = run()

# function 
# end

function loss(params)
    _prob = remake(prob, p=params);
    sol = solve(_prob, RK4(), tstops=tstops, saveat=tstops, adaptive=false, sensealg=BacksolveAdjoint());
    return -sum([sol.u[i].u.x[3][2] for i in 1:length(sol.t)])
end

params = [0.05; 0.05];
dp = Zygote.gradient(loss, params);
