using GLMakie
using Observables
using DifferentialEquations
using LinearAlgebra

# ----------------------------- System Parameters -----------------------------
const mpivot = 15.0    # mass of the pivot/rocket
const mp = 5.0         # mass of the pendulum
const l = 1.0          # length of the pendulum rod
const g = 9.81         # gravitational acceleration
const Bpivot = 5.0     # pivot damping coefficient
const Bθ = 0.5         # pendulum azimuthal damping coefficient
const Bφ = 0.5         # pendulum polar damping coefficient
const mass_decay_rate = 0.01  # rate at which rocket mass decreases (kg/s) - reduced for stability

struct Parameters
    mpivot::Float64  
    mp::Float64  
    l::Float64  
    g::Float64   
    Bpivot::Float64  
    Bθ::Float64
    Bφ::Float64
    mass_decay_rate::Float64
end

params = Parameters(mpivot, mp, l, g, Bpivot, Bθ, Bφ, mass_decay_rate)

# ----------------------------- 3D Rocket-Pendulum Model -----------------------------
function rocket_pendulum_3d!(dz, z, p, t)
    # Unpack state variables
    # z = [x, y, z, θ, φ, ẋ, ẏ, ż, θ̇, φ̇, mpivot]
    x, y, z_pos, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot, m_pivot = z
    
    # Enforce constraints to prevent numerical instability
    # Limit velocities and rates to reasonable values
    x_dot = clamp(x_dot, -10.0, 10.0)
    y_dot = clamp(y_dot, -10.0, 10.0)
    z_dot = clamp(z_dot, -10.0, 10.0)
    θ_dot = clamp(θ_dot, -5.0, 5.0)
    φ_dot = clamp(φ_dot, -5.0, 5.0)
    
    # Ensure mass stays positive and reasonable
    m_pivot = max(m_pivot, 0.1)
    
    # Calculate time-varying rocket mass
    m_pivot_dot = -p.mass_decay_rate  # mass decreases at a constant rate
    
    # Stop mass reduction if mass is too low
    if m_pivot <= 0.5
        m_pivot_dot = 0.0
    end
    
    # Control inputs (can be replaced with actual control law)
    Fx = 0.0  # External force in x
    Fy = 0.0  # External force in y
    Fz = m_pivot * p.g + 2.0 * sin(0.5 * t)  # Basic thrust to counteract gravity + smaller oscillation
    
    # Nonlinear disturbance (example: wind force proportional to velocity squared)
    # Reduced magnitude for stability
    Dx = 0.05 * x_dot^2 * (x_dot > 0 ? 1.0 : -1.0)
    Dy = 0.05 * y_dot^2 * (y_dot > 0 ? 1.0 : -1.0)
    Dz = 0.05 * z_dot^2 * (z_dot > 0 ? 1.0 : -1.0)
    
    # Calculate pendulum position relative to pivot
    x_rel = p.l * sin(φ) * cos(θ)
    y_rel = p.l * sin(φ) * sin(θ)
    z_rel = p.l * cos(φ)
    
    # Calculate pendulum acceleration terms
    # Terms for θ̈ with safety checks for numerical stability
    θ_accel_term = 0.0
    if abs(sin(φ)) > 1e-3  # Use a larger threshold for better stability
        θ_accel_term = p.l * sin(φ) * (
            -2.0 * φ_dot * θ_dot * cos(φ) / sin(φ) - 
            (x_dot * sin(θ) - y_dot * cos(θ)) / (p.l * sin(φ))
        )
    end
    
    # Terms for φ̈
    φ_accel_term = -p.l * (
        φ_dot^2 * sin(φ) * cos(φ) + 
        θ_dot^2 * sin(φ) * cos(φ) + 
        (x_dot * cos(θ) + y_dot * sin(θ)) / p.l
    )
    
    # Build the mass matrix for the coupled system
    # Order: [ẍ, ÿ, z̈, θ̈, φ̈]
    M = zeros(5, 5)
    
    # Add small regularization terms to the diagonal for numerical stability
    ε = 1e-6
    
    # Mass matrix components for rocket position
    M[1,1] = m_pivot + p.mp + ε
    M[2,2] = m_pivot + p.mp + ε
    M[3,3] = m_pivot + p.mp + ε
    
    # Cross-coupling between rocket and pendulum
    M[1,4] = -p.mp * p.l * sin(φ) * sin(θ)
    M[1,5] = p.mp * p.l * cos(φ) * cos(θ)
    M[2,4] = p.mp * p.l * sin(φ) * cos(θ)
    M[2,5] = p.mp * p.l * cos(φ) * sin(θ)
    M[3,5] = -p.mp * p.l * sin(φ)
    
    # Mass matrix components for pendulum angles
    M[4,1] = -p.mp * p.l * sin(φ) * sin(θ)
    M[4,2] = p.mp * p.l * sin(φ) * cos(θ)
    M[4,4] = p.mp * p.l^2 * sin(φ)^2 + ε
    
    M[5,1] = p.mp * p.l * cos(φ) * cos(θ)
    M[5,2] = p.mp * p.l * cos(φ) * sin(θ)
    M[5,3] = -p.mp * p.l * sin(φ)
    M[5,5] = p.mp * p.l^2 + ε
    
    # Right-hand side vector
    B = zeros(5)
    
    # Forces on rocket (including gravity, control, and disturbances)
    B[1] = Fx - Dx - p.Bpivot * x_dot - m_pivot_dot * x_dot
    B[2] = Fy - Dy - p.Bpivot * y_dot - m_pivot_dot * y_dot
    B[3] = Fz - Dz - p.Bpivot * z_dot - (m_pivot + p.mp) * p.g - m_pivot_dot * z_dot
    
    # Pendulum dynamics (including gravity and damping)
    B[4] = p.mp * p.g * p.l * sin(φ) * sin(θ) - p.Bθ * θ_dot + θ_accel_term
    B[5] = p.mp * p.g * p.l * sin(φ) * cos(θ) - p.Bφ * φ_dot + φ_accel_term
    
    # Solve the system using linear algebra
    accelerations = M \ B
    
    # Apply stability constraints to accelerations
    accelerations = clamp.(accelerations, -20.0, 20.0)
    
    # Set the derivatives
    dz[1] = x_dot                # ẋ
    dz[2] = y_dot                # ẏ
    dz[3] = z_dot                # ż
    dz[4] = θ_dot                # θ̇
    dz[5] = φ_dot                # φ̇
    dz[6] = accelerations[1]     # ẍ
    dz[7] = accelerations[2]     # ÿ
    dz[8] = accelerations[3]     # z̈
    dz[9] = accelerations[4]     # θ̈
    dz[10] = accelerations[5]    # φ̈
    dz[11] = m_pivot_dot         # ṁpivot
end

# ----------------------- Create the ODE Problem & Solve ----------------------
# Initial conditions: [x, y, z, θ, φ, ẋ, ẏ, ż, θ̇, φ̇, mpivot]
z0 = [0.0, 0.0, 0.0, 0.1, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, mpivot] 
tspan = (0.0, 35.0)

# Solve the model ODEs with a more robust solver for stiff systems
prob = ODEProblem(rocket_pendulum_3d!, z0, tspan, params)
sol = solve(prob, Tsit5(), abstol=1e-6, reltol=1e-6, maxiters=1000000)

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------- Create Figure + Plots -----------------------------
# Create 3D visualization
fig1 = Figure(size=(800, 800), fontsize=12)
ax = fig1[1, 1] = Axis3(fig1, 
                     xlabel = "x", ylabel = "y", zlabel = "z",
                     limits = (-2*l, 2*l, -2*l, 2*l, -2*l, 2*l),
                     aspect = :data)

# Initial positions
rocket_pos = Point3f(z0[1], z0[2], z0[3])
pendulum_pos = Point3f(
    z0[1] + l * sin(z0[5]) * cos(z0[4]), 
    z0[2] + l * sin(z0[5]) * sin(z0[4]), 
    z0[3] + l * cos(z0[5])
)

# Create visualization elements (using meshscatter for better compatibility)
rocket_plot = meshscatter!(ax, [rocket_pos[1]], [rocket_pos[2]], [rocket_pos[3]], markersize = 0.2, color = :red)
rod_plot = lines!(ax, [rocket_pos[1], pendulum_pos[1]], [rocket_pos[2], pendulum_pos[2]], [rocket_pos[3], pendulum_pos[3]], linewidth = 3, color = :blue)
pendulum_plot = meshscatter!(ax, [pendulum_pos[1]], [pendulum_pos[2]], [pendulum_pos[3]], markersize = 0.15, color = :blue)

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
                      title = "Rocket Height")

# Create empty line objects for trajectories
rocket_traj = lines!(ax_3d_traj, Float64[], Float64[], Float64[], color = :red)
pendulum_traj = lines!(ax_3d_traj, Float64[], Float64[], Float64[], color = :blue)

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
fps = 60
dt_frame = 1/fps
t_end = sol.t[end]

@async begin
    t_sim = 0.0
    
    while t_sim <= t_end && t_sim <= sol.t[end]
        try
            # Sample the solution at the current simulation timestep
            z = sol(t_sim)
            
            # Unpack state
            x, y, z_height_val, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot, m_pivot = z
            
            # Calculate pendulum position
            x_pend = x + l * sin(φ) * cos(θ)
            y_pend = y + l * sin(φ) * sin(θ)
            z_pend = z_height_val + l * cos(φ)
            
            # Update rocket size based on remaining mass
            rocket_size = 0.1 + 0.1 * (m_pivot / mpivot)
            
            # Update visualization elements
            rocket_plot[1] = [x]
            rocket_plot[2] = [y]
            rocket_plot[3] = [z_height_val]
            rocket_plot.markersize = rocket_size
            
            rod_plot[1] = [x, x_pend]
            rod_plot[2] = [y, y_pend]
            rod_plot[3] = [z_height_val, z_pend]
            
            pendulum_plot[1] = [x_pend]
            pendulum_plot[2] = [y_pend]
            pendulum_plot[3] = [z_pend]
            
            # Store trajectory data
            push!(rocket_x, x)
            push!(rocket_y, y)
            push!(rocket_z, z_height_val)
            push!(pendulum_x, x_pend)
            push!(pendulum_y, y_pend)
            push!(pendulum_z, z_pend)
            
            # Store phase portrait data
            push!(theta, θ)
            push!(thetadot, θ_dot)
            push!(phi, φ)
            push!(phidot, φ_dot)
            push!(time_array, t_sim)
            push!(z_height, z_height_val)
            
            # Update trajectory visualizations
            rocket_traj[1] = rocket_x
            rocket_traj[2] = rocket_y
            rocket_traj[3] = rocket_z
            pendulum_traj[1] = pendulum_x
            pendulum_traj[2] = pendulum_y
            pendulum_traj[3] = pendulum_z
            
            # Update phase portraits
            theta_line[1] = theta
            theta_line[2] = thetadot
            phi_line[1] = phi
            phi_line[2] = phidot
            z_line[1] = time_array
            z_line[2] = z_height
            
            # Update Observables
            time_obs[] = t_sim
            state_obs[] = z
            
            sleep(dt_frame)
            t_sim += dt_frame
        catch e
            println("Error at t=$t_sim: $e")
            break
        end
    end
end

println("3D Rocket-Pendulum simulation is running!")