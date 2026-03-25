using DifferentialEquations, Plots

# --- Define the Model Function ---
function consumer_resource!(du, u, p, t)
    # Unpack parameters and dimensions
    L_r, L_s, V_max, K_s, Y, m, D, periods = p
    
    # Split state vector: Resources first, then Species
    R = @view u[1:L_r]
    N = @view u[L_r+1:end]
    
    dR = @view du[1:L_r]
    dN = @view du[L_r+1:end]

    # 1. Calculate per-species growth rates across all resources
    # growth_matrix[i, j] = V_max_ij * R_j / (K_s_ij + R_j)
    growth_matrix = V_max .* (R' ./ (K_s .+ R'))
    
    # Total growth for each species i (sum across resources j)
    total_growth = sum(growth_matrix, dims=2)[:] 

    # 2. Species Dynamics (dN/dt)
    for i in 1:L_s
        dN[i] = N[i] * (total_growth[i] - m[i])
    end

    # 3. Resource Dynamics (dR/dt)
    for j in 1:L_r
        # Cyclic supply for resource j
        S_in = 0.5 + 0.4 * sin(2π * t / periods[j])
        
        # Total consumption of resource j by all species i
        consumption_j = sum((1 ./ Y[:, j]) .* growth_matrix[:, j] .* N)
        
        dR[j] = D * (S_in - R[j]) - consumption_j
    end
end

# --- Setup for ANY L_s and L_r ---
L_r = 4  # Number of Resources
L_s = 6  # Number of Species

# Generate random but biologically plausible parameters
p = (
    L_r = L_r,
    L_s = L_s,
    V_max = rand(L_s, L_r) .* 0.5,     # Max growth rates
    K_s   = rand(L_s, L_r) .* 0.2,     # Half-saturation
    Y     = rand(L_s, L_r) .* 0.4 .+ 0.2, # Yields
    m     = rand(L_s) .* 0.1,          # Mortalities
    D     = 0.1,                       # Dilution
    periods = rand(L_r) .* 20 .+ 10    # Unique supply cycle per resource
)

# Initial conditions: random concentrations
u0 = [rand(L_r); rand(L_s)]
tspan = (0.0, 300.0)

# Solve
prob = ODEProblem(consumer_resource!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-6)

# --- Plot Results ---
p1 = plot(sol, vars=1:L_r, title="Resources ($L_r)", ylabel="R", palette=:blues)
p2 = plot(sol, vars=L_r+1:L_r+L_s, title="Species ($L_s)", ylabel="N", palette=:greens)
plot(p1, p2, layout=(2,1), size=(900, 700), legend=false)
