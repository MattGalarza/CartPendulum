using GLMakie, Observables, DifferentialEquations, LinearAlgebra, Sundials, OrdinaryDiffEq

# ----------------------------- System Parameters -----------------------------
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

# Zero-control and zero-disturbance functions for now
zero_ctrl(t, u) = zeros(3)
zero_disturbance(t, u) = zeros(3)

# Instantiate parameters
params = Parameters(
    5.0,    # mpivot
    5.0,    # mp
    1.0,    # ℓ
    9.81,   # g
    1.0,    # Bpivot
    1.0,    # Bp
    zero_ctrl,
    zero_disturbance
)

# ----------------------------- Quaternion Utilities -----------------------------
# Convert Euler angles (θ, φ) to quaternion [qw, qx, qy, qz]
function euler_to_quaternion(θ, φ)
    # For 3D pendulum: φ is polar angle from z, θ is azimuthal angle in xy-plane
    # Rotation around z-axis by θ, then around the new y-axis by φ
    
    # Calculate half-angles
    θ_2 = θ/2
    φ_2 = φ/2
    
    # Calculate quaternion components
    qw = cos(φ_2) * cos(θ_2)
    qx = sin(φ_2) * cos(θ_2)
    qy = sin(φ_2) * sin(θ_2)
    qz = cos(φ_2) * sin(θ_2)
    
    return [qw, qx, qy, qz]
end

# Convert quaternion to direction vector (unit vector along pendulum)
function quaternion_to_direction(q)
    qw, qx, qy, qz = q
    
    # Transform unit z-vector using quaternion rotation
    # This gives us the direction the pendulum is pointing
    dir_x = 2.0 * (qx*qz + qw*qy)
    dir_y = 2.0 * (qy*qz - qw*qx)
    dir_z = 1.0 - 2.0 * (qx*qx + qy*qy)
    
    return [dir_x, dir_y, dir_z]
end

# Function to normalize a quaternion
function normalize_quaternion(q)
    norm = sqrt(sum(q.^2))
    return q ./ norm
end

# ----------------------------- Dynamics Function -----------------------------
function pendulum_quaternion!(du, u, p, t)
    # Unpack state: [x, y, z, qw, qx, qy, qz, x_dot, y_dot, z_dot, ω_x, ω_y, ω_z]
    # where q is quaternion, ω is angular velocity in body frame
    x, y, z = u[1:3]
    q = u[4:7]
    x_dot, y_dot, z_dot = u[8:10]
    ω = u[11:13]  # angular velocity in body frame
    
    # Normalize quaternion to prevent drift
    q = normalize_quaternion(q)
    
    # Get the direction of the pendulum (unit vector)
    dir = quaternion_to_direction(q)
    
    # Calculate pendulum endpoint position
    xp = x + p.l * dir[1]
    yp = y + p.l * dir[2]
    zp = z + p.l * dir[3]
    
    # Calculate pendulum velocity (analytical from position and angular velocity)
    # Cross product of angular velocity with pendulum vector gives tangential velocity
    ω_cross_dir = [
        ω[2]*dir[3] - ω[3]*dir[2],
        ω[3]*dir[1] - ω[1]*dir[3],
        ω[1]*dir[2] - ω[2]*dir[1]
    ]
    
    xp_dot = x_dot + p.l * ω_cross_dir[1]
    yp_dot = y_dot + p.l * ω_cross_dir[2]
    zp_dot = z_dot + p.l * ω_cross_dir[3]
    
    # External inputs (Control and Disturbance)
    # Adapt state vector for compatibility with control functions
    u_adapted = [x, y, z, 0, 0, x_dot, y_dot, z_dot, 0, 0]  # Placeholder for θ, φ
    F = p.control(t, u_adapted)
    Dp = p.disturbance(t, u_adapted)
    
    # Drag (nonlinear damping ∝ v^2)
    drag_pend = p.Bp * [
        sign(xp_dot) * xp_dot^2,
        sign(yp_dot) * yp_dot^2,
        sign(zp_dot) * zp_dot^2
    ]
    
    drag_pivot = p.Bpivot * [
        sign(x_dot) * x_dot^2,
        sign(y_dot) * y_dot^2,
        sign(z_dot) * z_dot^2
    ]
    
    # Build the translational equations (for pivot)
    # F = ma + drag + reaction_force
    m_total = p.mpivot + p.mp
    
    # Calculate reaction force from pendulum (through rod tension)
    # This involves the centripetal acceleration and gravity components
    centripetal_accel = [
        p.l * (ω[2]^2 + ω[3]^2) * dir[1] - p.l * ω[1] * (ω[2]*dir[2] + ω[3]*dir[3]),
        p.l * (ω[1]^2 + ω[3]^2) * dir[2] - p.l * ω[2] * (ω[1]*dir[1] + ω[3]*dir[3]),
        p.l * (ω[1]^2 + ω[2]^2) * dir[3] - p.l * ω[3] * (ω[1]*dir[1] + ω[2]*dir[2])
    ]
    
    # Solve for translational accelerations
    pivot_accel = [
        (F[1] - drag_pivot[1] - drag_pend[1] + p.mp * centripetal_accel[1]) / m_total,
        (F[2] - drag_pivot[2] - drag_pend[2] + p.mp * centripetal_accel[2]) / m_total,
        (F[3] - drag_pivot[3] - drag_pend[3] + p.mp * centripetal_accel[3] - m_total * p.g) / m_total
    ]
    
    # Calculate torque on pendulum (cross product of position vector and force)
    gravity_force = [0, 0, -p.mp * p.g]
    r_vec = p.l .* dir  # position vector from pivot to pendulum mass
    
    # Torque = r × F for gravity and drag forces
    gravity_torque = [
        r_vec[2]*gravity_force[3] - r_vec[3]*gravity_force[2],
        r_vec[3]*gravity_force[1] - r_vec[1]*gravity_force[3],
        r_vec[1]*gravity_force[2] - r_vec[2]*gravity_force[1]
    ]
    
    drag_torque = [
        r_vec[2]*drag_pend[3] - r_vec[3]*drag_pend[2],
        r_vec[3]*drag_pend[1] - r_vec[1]*drag_pend[3],
        r_vec[1]*drag_pend[2] - r_vec[2]*drag_pend[1]
    ]
    
    # Calculate angular acceleration from torque
    # For a simple pendulum with mass at the end of a rod, inertia tensor is:
    # I = mp * l² * (identity_matrix - direction_outer_product)
    inertia_tensor = p.mp * p.l^2 * (Matrix{Float64}(I, 3, 3) - dir * dir')
    
    # Add a small regularization to ensure invertibility
    inertia_tensor += 1e-6 * Matrix{Float64}(I, 3, 3)
    
    # Calculate angular acceleration: α = I⁻¹ * (τ - ω × (I·ω))
    angular_momentum = inertia_tensor * ω
    gyroscopic_torque = [
        ω[2]*angular_momentum[3] - ω[3]*angular_momentum[2],
        ω[3]*angular_momentum[1] - ω[1]*angular_momentum[3],
        ω[1]*angular_momentum[2] - ω[2]*angular_momentum[1]
    ]
    
    total_torque = gravity_torque - drag_torque - gyroscopic_torque
    angular_accel = inertia_tensor \ total_torque
    
    # Quaternion derivative from angular velocity
    q_dot = 0.5 * [
        -q[2]*ω[1] - q[3]*ω[2] - q[4]*ω[3],
         q[1]*ω[1] + q[3]*ω[3] - q[4]*ω[2],
         q[1]*ω[2] - q[2]*ω[3] + q[4]*ω[1],
         q[1]*ω[3] + q[2]*ω[2] - q[3]*ω[1]
    ]
    
    # Update the state derivative vector
    # Position derivatives (velocities)
    du[1:3] = [x_dot, y_dot, z_dot]
    
    # Quaternion derivatives
    du[4:7] = q_dot
    
    # Velocity derivatives (accelerations)
    du[8:10] = pivot_accel
    
    # Angular velocity derivatives (angular accelerations)
    du[11:13] = angular_accel
end

# ----------------------------- Set up & Solve ODE Problem -----------------------------
# Convert initial spherical coordinates to quaternion
θ_init = 0.1
φ_init = 0.3
q_init = euler_to_quaternion(θ_init, φ_init)

# Initial condition: [x, y, z, qw, qx, qy, qz, x_dot, y_dot, z_dot, ω_x, ω_y, ω_z]
z0 = [0.0, 0.0, 0.0, q_init..., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
tspan = (0.0, 35.0)

# Solve the pendulum system
prob = ODEProblem(pendulum_quaternion!, z0, tspan, params)
sol = solve(prob, CVODE_BDF(), abstol=1e-8, reltol=1e-8, maxiters=1000000)

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------- Post-processing -----------------------------
# Extract positions for visualization
num_points = min(1000, length(sol.t))
indices = round.(Int, range(1, length(sol.t), length=num_points))

# Arrays to store positions
pivot_positions = zeros(num_points, 3)
pendulum_positions = zeros(num_points, 3)

for (i, idx) in enumerate(indices)
    u = sol.u[idx]
    x, y, z = u[1:3]
    q = normalize_quaternion(u[4:7])
    dir = quaternion_to_direction(q)
    
    pivot_positions[i, :] = [x, y, z]
    pendulum_positions[i, :] = [x + params.l * dir[1], 
                               y + params.l * dir[2], 
                               z + params.l * dir[3]]
end

# Now you can use these positions for visualization
println("Solution successfully processed for visualization")

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
            u = sol(t_sim)
            
            # Unpack state
            x, y, z, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot = u
            
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
