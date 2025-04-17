using GLMakie, Observables, DifferentialEquations, LinearAlgebra, Sundials, OrdinaryDiffEq

# ----------------------------- System Parameters -----------------------------
# Now includes two Function fields: one for control forces, one for disturbances
struct Parameters
    mpivot::Float64    # pivot (rocket) mass
    mp::Float64        # pendulum mass
    l::Float64         # rod length
    g::Float64         # gravitational acceleration
    Bpivot::Float64    # pivot drag coefficient (½ρC_dA)
    Bp::Float64        # pendulum drag coefficient (½ρC_dA)
    control::Function  # (t,u) -> [Fx, Fy, Fz]
    disturbance::Function  # (t,u) -> [Dp_x, Dp_y, Dp_z]
end

# Zero‐control and zero‐disturbance functions for now
zero_ctrl(t, u) = zeros(3)
zero_disturbance(t, u) = zeros(3)

# Instantiate parameters with zeros for now
params = Parameters(
    5.0,    # mpivot
    5.0,    # mp
    1.0,    # ℓ
    9.81,   # g
    5.0,    # Bpivot
    2.0,    # Bp
    zero_ctrl,
    zero_disturbance
)

# ----------------------------- Dynamics Function -----------------------------
function pendulum_3d!(du, u, p, t)
    # Unpack state: [x, y, z, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot]
    x, y, z, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot = u

    # Kinematics of pendulum bob
    xp  = x  + (p.l * sin(φ) * cos(θ))
    yp  = y  + (p.l * sin(φ) * sin(θ))
    zp  = z  + (p.l * cos(φ))
    xp_dot = x_dot + (p.l * (φ_dot * cos(φ) * cos(θ) - θ_dot * sin(φ) * sin(θ)))
    yp_dot = y_dot + (p.l * (φ_dot * cos(φ) * sin(θ) + θ_dot * sin(φ) * cos(θ)))
    zp_dot = z_dot - (p.l * (φ_dot * sin(φ)))

    # External inputs (Contol and Disturbance)
    F   = p.control(t, u)        # [Fx, Fy, Fz]
    Dp  = p.disturbance(t, u)    # [Dp_x, Dp_y, Dp_z]

    # Drag (nonlinear damping ∝ v^2)
    drag_pend  = p.Bp * [xp_dot^2, yp_dot^2, zp_dot^2]
    drag_pivot = p.Bpivot * [x_dot^2,  y_dot^2,  z_dot^2]

    # Build mass matrix M and RHS vector B
    T = eltype(u)
    M = zeros(T, 5, 5)
    B = zeros(T, 5)

    # Add a small regularization to prevent M matrix from becoming singular
    ϵ = 1e-6

    # Pivot + pendulum diagonal (ẍ, ÿ, z̈)
    M[1,1] = p.mpivot + p.mp + ϵ
    M[2,2] = p.mpivot + p.mp + ϵ
    M[3,3] = p.mpivot + p.mp + ϵ

    # Geometric coupling entries
    M[1,4] = -p.mp * p.l*sin(φ)*sin(θ)
    M[1,5] =  p.mp * p.l*cos(φ)*cos(θ)
    M[2,4] =  p.mp * p.l*sin(φ)*cos(θ)
    M[2,5] =  p.mp * p.l*cos(φ)*sin(θ)
    M[3,5] = -p.mp * p.l*sin(φ)
    # symmetry
    M[4,1] = M[1,4]
    M[5,1] = M[1,5]
    M[4,2] = M[2,4]
    M[5,2] = M[2,5]
    M[5,3] = M[3,5]

    # Pendulum angular inertias
    M[4,4] = p.mp * p.l^2 * sin(φ)^2 + ϵ
    M[5,5] = p.mp * p.l^2 + ϵ

    # Fill force vector B
    # Pivot translation equations
    B[1] = F[1] - drag_pivot[1] - Dp[1]
    B[2] = F[2] - drag_pivot[2] - Dp[2]
    B[3] = F[3] - drag_pivot[3] - (p.mpivot + p.mp)*p.g

    # Subtract pendulum drag & disturbance (reaction on pivot)
    B[1] += -drag_pend[1] - Dp[1]
    B[2] += -drag_pend[2] - Dp[2]
    B[3] += -drag_pend[3] - Dp[3]

    # Pendulum angular “forcing” (gravity torques)
    B[4] = -p.mp * p.g * p.l * sin(φ)
    B[5] = -p.mp * p.g * p.l * sin(φ) * cos(θ)

    # Solve for accelerations
    accels = M \ B

    # hack: limit accelerations to ±Amax
    #Amax = 500.0
    #accels = clamp.(accels, -Amax, Amax)

    # Write the derivatives
    du[1]  = x_dot
    du[2]  = y_dot
    du[3]  = z_dot
    du[4]  = θ_dot
    du[5]  = φ_dot
    du[6]  = accels[1]   # ẍ
    du[7]  = accels[2]   # ÿ
    du[8]  = accels[3]   # z̈
    du[9]  = accels[4]   # θ̈
    du[10] = accels[5]   # φ̈
end

# ----------------------- Set up & Solve ODE Problem -------------------------
# Initial condition: [x, y, z, θ, φ, ẋ, ẏ, ż, θ̇, φ̇]
z0    = [0.0, 0.0, 0.0, 0.1, 0.3, 0, 0, 0, 0, 0]
tspan = (0.0, 35.0)

# Solve the pendulum system
prob = ODEProblem(pendulum_3d!, z0, tspan, params)
sol  = solve(prob, Rodas5(), abstol=1e-8, reltol=1e-10, maxiters=1000000, force_dtmin=true)

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------- Create Figure + Plots -----------------------------
# Create 3D visualization
fig1 = Figure(size=(800, 800), fontsize=12)
ax = fig1[1, 1] = Axis3(fig1, 
                     xlabel = "x", ylabel = "y", zlabel = "z",
                     limits = (-5*params.l, 5*params.l, -5*params.l, 5*params.l, -5*params.l, 5*params.l),
                     aspect = :data)

# Initial positions
rocket_pos = Point3f(z0[1], z0[2], z0[3])
pendulum_pos = Point3f(
    z0[1] + params.l * sin(z0[5]) * cos(z0[4]), 
    z0[2] + params.l * sin(z0[5]) * sin(z0[4]), 
    z0[3] + params.l * cos(z0[5])
)

# Create visualization elements
rocket_plot = meshscatter!(ax, [rocket_pos[1]], [rocket_pos[2]], [rocket_pos[3]], 
                          markersize = 0.2, color = :red)
rod_plot = lines!(ax, [rocket_pos[1], pendulum_pos[1]], [rocket_pos[2], pendulum_pos[2]], 
                 [rocket_pos[3], pendulum_pos[3]], linewidth = 3, color = :blue)
pendulum_plot = meshscatter!(ax, [pendulum_pos[1]], [pendulum_pos[2]], [pendulum_pos[3]], 
                            markersize = 0.15, color = :blue)

display(GLMakie.Screen(), fig1)

# Create trajectory visualization
fig2 = Figure(size=(1200, 800))
ax_3d_traj = fig2[1, 1] = Axis3(fig2, 
                            xlabel = "x", ylabel = "y", zlabel = "z",
                            title = "3D Trajectory",
                            aspect = :data)
                            
# Create phase portrait plots
ax_theta = fig2[1, 2] = Axis(fig2, 
                          xlabel = "θ (rad)", ylabel = "θ̇ (rad/s)", 
                          title = "Azimuthal Phase Portrait")
ax_phi = fig2[2, 1] = Axis(fig2, 
                        xlabel = "φ (rad)", ylabel = "φ̇ (rad/s)", 
                        title = "Polar Phase Portrait")
ax_z = fig2[2, 2] = Axis(fig2, 
                      xlabel = "t (s)", ylabel = "z (m)", 
                      title = "Pivot Height")

# Create empty line objects for trajectories
rocket_traj = lines!(ax_3d_traj, Float64[], Float64[], Float64[], color = :red, label = "Pivot")
pendulum_traj = lines!(ax_3d_traj, Float64[], Float64[], Float64[], color = :blue, label = "Pendulum")
Legend(fig2[1, 3], ax_3d_traj)

theta_line = lines!(ax_theta, Float64[], Float64[], color = :purple)
phi_line = lines!(ax_phi, Float64[], Float64[], color = :green)
z_line = lines!(ax_z, Float64[], Float64[], color = :orange)

display(GLMakie.Screen(), fig2)

# Observable state for animation
state_obs = Observable(z0)
time_obs = Observable(0.0)

# Storage for trajectories and phase diagrams
rocket_x = Float64[]
rocket_y = Float64[]
rocket_z = Float64[]
pendulum_x = Float64[]
pendulum_y = Float64[]
pendulum_z = Float64[]
theta = Float64[]
thetadot = Float64[]
phi = Float64[]
phidot = Float64[]
time_array = Float64[]
z_height = Float64[]

# ----------------------------- Update Function + visual -----------------------------
fps = 120
dt_frame = 1/fps
t_end = sol.t[end]

# Initialize trajectory arrays (define them here to ensure correct scope)
rocket_x = Float64[]
rocket_y = Float64[]
rocket_z = Float64[]
pendulum_x = Float64[]
pendulum_y = Float64[]
pendulum_z = Float64[]
theta = Float64[]
thetadot = Float64[]
phi = Float64[]
phidot = Float64[]
time_array = Float64[]
z_height = Float64[]

# Create animation focusing only on pendulum movement, no trails
sleep(3.0) # Gives delay to visualization
@async begin
    t_sim = 0.0
    
    while t_sim <= t_end && t_sim <= sol.t[end]
        try
            # Sample the solution at the current simulation timestep
            z = sol(t_sim)
            
            # Unpack state
            x, y, z_height_val, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot = z
            
            # Calculate pendulum position
            x_pend = x + params.l * sin(φ) * cos(θ)
            y_pend = y + params.l * sin(φ) * sin(θ)
            z_pend = z_height_val + params.l * cos(φ)
            
            # Update visualization elements - show actual movement
            rocket_plot[1] = [x]
            rocket_plot[2] = [y]
            rocket_plot[3] = [z_height_val]
            
            rod_plot[1] = [x, x_pend]
            rod_plot[2] = [y, y_pend]
            rod_plot[3] = [z_height_val, z_pend]
            
            pendulum_plot[1] = [x_pend]
            pendulum_plot[2] = [y_pend]
            pendulum_plot[3] = [z_pend]
            
            # Update Observables
            time_obs[] = t_sim
            state_obs[] = z
            
            sleep(dt_frame)
            t_sim += dt_frame
        catch e
            println("Error at t=$t_sim: $e")
            println("Error type: ", typeof(e))
            break
        end
    end
end

println("3D Pendulum simulation is running!")



