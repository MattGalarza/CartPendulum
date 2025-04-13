using GLMakie
using Observables
using DifferentialEquations
using LinearAlgebra

# ----------------------------- System Parameters -----------------------------
const mpivot = 15.0    # mass of the pivot/rocket
const mp = 5.0         # mass of the pendulum
const l = 1.0          # length of the pendulum rod
const g = 9.81         # gravitational acceleration
const Bpivot = 2.0     # pivot damping coefficient
const Bp = 2.0         # pendulum damping coefficient
const disturbance_amplitude = 0.0  # amplitude of disturbance forces

struct Parameters
    mpivot::Float64  
    mp::Float64  
    l::Float64  
    g::Float64   
    Bpivot::Float64  
    Bp::Float64
    disturbance_amplitude::Float64
end

params = Parameters(mpivot, mp, l, g, Bpivot, Bp, disturbance_amplitude)

# ----------------------------- 3D Pendulum Model -----------------------------
function pendulum_3d!(dz, z, p, t)
    # Unpack state variables
    # z = [x, y, z, θ, φ, ẋ, ẏ, ż, θ̇, φ̇]
    x, y, z_pos, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot = z
    
    # Calculate pendulum position and velocity relative to pivot
    # Position
    x_p = x + p.l * sin(φ) * cos(θ)
    y_p = y + p.l * sin(φ) * sin(θ)
    z_p = z_pos + p.l * cos(φ)
    
    # Velocity
    x_p_dot = x_dot + p.l * (θ_dot * cos(φ) * cos(θ) - φ_dot * sin(φ) * sin(θ))
    y_p_dot = y_dot + p.l * (θ_dot * cos(φ) * sin(θ) + φ_dot * sin(φ) * cos(θ))
    z_p_dot = z_dot - p.l * φ_dot * sin(φ)
    
    # Generate disturbances (time-varying sinusoidal disturbances)
    Dp_x = p.disturbance_amplitude * sin(0.5*t)
    Dp_y = p.disturbance_amplitude * cos(0.7*t)
    Dp_z = p.disturbance_amplitude * sin(0.3*t)
    
    # Control forces (can be replaced with actual control law)
    Fx = 0.0  # No control for now
    Fy = 0.0  # No control for now
    Fz = p.mpivot * p.g  # Just enough to counter gravity
    
    # Build the mass matrix for the coupled system
    # Order: [ẍ, ÿ, z̈, θ̈, φ̈]
    M = zeros(5, 5)
    
    # Add small regularization terms to the diagonal for numerical stability
    ε = 1e-6
    
    # Mass matrix components for pivot position
    M[1,1] = p.mpivot + p.mp + ε
    M[2,2] = p.mpivot + p.mp + ε
    M[3,3] = p.mpivot + p.mp + ε
    
    # Cross-coupling between pivot and pendulum
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
    
    # Pendulum damping terms
    B_pend_x = p.Bp * x_p_dot
    B_pend_y = p.Bp * y_p_dot
    B_pend_z = p.Bp * z_p_dot
    
    # Pivot damping terms
    B_pivot_x = p.Bpivot * x_dot
    B_pivot_y = p.Bpivot * y_dot
    B_pivot_z = p.Bpivot * z_dot
    
    # Right-hand side vector (B) with all forces
    B = zeros(5)
    
    # Terms in acceleration equations
    θ_accel_term = p.l * sin(φ) * (
        -2.0 * φ_dot * θ_dot * cos(φ) / max(sin(φ), 1e-3) - 
        (x_dot * sin(θ) - y_dot * cos(θ)) / (p.l * max(sin(φ), 1e-3))
    )
    
    φ_accel_term = -p.l * (
        φ_dot^2 * sin(φ) * cos(φ) + 
        θ_dot^2 * sin(φ) * cos(φ) + 
        (x_dot * cos(θ) + y_dot * sin(θ)) / p.l
    )
    
    # Forces affecting pivot accelerations
    B[1] = Fx - B_pivot_x + Dp_x
    B[2] = Fy - B_pivot_y + Dp_y
    B[3] = Fz - B_pivot_z - (p.mpivot + p.mp) * p.g
    
    # Add pendulum reaction forces (negative of pendulum acceleration terms)
    B[1] += p.mp * (
        -p.l * (θ̈_term = θ_dot^2 * sin(φ) * cos(θ) + 2 * θ_dot * φ_dot * cos(φ) * sin(θ) + φ_dot^2 * sin(φ) * cos(θ))
    ) + B_pend_x
    
    B[2] += p.mp * (
        -p.l * (θ̈_term = θ_dot^2 * sin(φ) * sin(θ) - 2 * θ_dot * φ_dot * cos(φ) * cos(θ) + φ_dot^2 * sin(φ) * sin(θ))
    ) + B_pend_y
    
    B[3] += p.mp * (
        p.l * (φ̈_term = φ_dot^2 * cos(φ))
    ) + B_pend_z
    
    # Forces affecting angular accelerations (including gravity and acceleration coupling)
    B[4] = p.mp * p.g * p.l * sin(φ) * sin(θ) + θ_accel_term - B_pend_x * p.l * sin(φ) * sin(θ) + B_pend_y * p.l * sin(φ) * cos(θ)
    B[5] = p.mp * p.g * p.l * sin(φ) * cos(θ) + φ_accel_term + B_pend_x * p.l * cos(φ) * cos(θ) + B_pend_y * p.l * cos(φ) * sin(θ) - B_pend_z * p.l * sin(φ)
    
    # Solve the system using linear algebra
    accelerations = M \ B
    
    # Apply stability constraints to accelerations to prevent numerical issues
    accelerations = clamp.(accelerations, -50.0, 50.0)
    
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
end

# ----------------------- Create the ODE Problem & Solve ----------------------
# Initial conditions: [x, y, z, θ, φ, ẋ, ẏ, ż, θ̇, φ̇]
# Starting with pendulum slightly off-center - not perfectly inverted
z0 = [0.0, 0.0, 0.0, 0.1, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0] 
tspan = (0.0, 35.0)

# Solve the model ODEs with higher resolution for smoother visualization
prob = ODEProblem(pendulum_3d!, z0, tspan, params)
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8, saveat=0.0001)

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------- Create Figure + Plots -----------------------------
# Create 3D visualization
fig1 = Figure(size=(800, 800), fontsize=12)
ax = fig1[1, 1] = Axis3(fig1, 
                     xlabel = "x", ylabel = "y", zlabel = "z",
                     limits = (-10*l, 10*l, -10*l, 10*l, -10*l, 10*l),
                     aspect = :data)

# Initial positions
rocket_pos = Point3f(z0[1], z0[2], z0[3])
pendulum_pos = Point3f(
    z0[1] + l * sin(z0[5]) * cos(z0[4]), 
    z0[2] + l * sin(z0[5]) * sin(z0[4]), 
    z0[3] + l * cos(z0[5])
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
@async begin
    t_sim = 0.0
    
    while t_sim <= t_end && t_sim <= sol.t[end]
        try
            # Sample the solution at the current simulation timestep
            z = sol(t_sim)
            
            # Unpack state
            x, y, z_height_val, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot = z
            
            # Calculate pendulum position
            x_pend = x + l * sin(φ) * cos(θ)
            y_pend = y + l * sin(φ) * sin(θ)
            z_pend = z_height_val + l * cos(φ)
            
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