using GLMakie, Observables, DifferentialEquations, LinearAlgebra, Sundials
using Rotations, StaticArrays, OrdinaryDiffEq

# ----------------------------- System Parameters -----------------------------
struct Parameters
    mpivot::Float64       # pivot (rocket) mass
    mp::Float64           # pendulum mass
    l::Float64            # rod length
    g::Float64            # gravity
    Bpivot::Float64       # pivot drag coefficient
    Bp::Float64           # pendulum drag coefficient
    control::Function     # (t,u) -> SVector{3}([F_x,F_y,F_z])
    disturbance::Function # (t,u) -> SVector{3}([Dp_x,Dp_y,Dp_z])
end

# zero‑input placeholders
zero_ctrl(t,u) = zero(SVector{3,eltype(u)})
zero_dist(t,u) = zero(SVector{3,eltype(u)})

params = Parameters(
    5.0,    # mpivot
    5.0,    # mp
    1.0,    # l
    9.81,   # g
    0.1,    # Bpivot
    0.1,    # Bp
    zero_ctrl,
    zero_dist
)

# ----------------------------- Dynamics Function -----------------------------
function pendulum_quat!(du, u, p::Parameters, t)
    # State: [ x, y, z,
    #          q0, q1, q2, q3,
    #          x_dot, y_dot, z_dot,
    #          omega_x, omega_y, omega_z ]
    x, y, z,
    q0, q1, q2, q3,
    x_dot, y_dot, z_dot,
    omega_x, omega_y, omega_z = u

    # Reconstruct quaternion and rod direction
    q = UnitQuaternion(q0, q1, q2, q3)
    n = q * Vec(0.0, 0.0, 1.0)    # unit‐vector along the rod

    # Bob position
    bob_x = x + p.l * n[1]
    bob_y = y + p.l * n[2]
    bob_z = z + p.l * n[3]

    # External inputs
    F   = p.control(t, u)       # SVector{3}
    Dp  = p.disturbance(t, u)   # SVector{3}

    # Bob velocity (translational + rotation)
    vp_x = x_dot + p.l*(omega_y*n[3] - omega_z*n[2])
    vp_y = y_dot + p.l*(omega_z*n[1] - omega_x*n[3])
    vp_z = z_dot + p.l*(omega_x*n[2] - omega_y*n[1])
    drag_bob = p.Bp * SVector(vp_x^2, vp_y^2, vp_z^2)

    # Pivot drag
    drag_pivot = p.Bpivot * SVector(x_dot^2, y_dot^2, z_dot^2)

    # -- Translational dynamics: M * accel = B --
    M_trans = (p.mpivot + p.mp) * I(3)
    B_trans = F .- drag_bob .- Dp .- SVector(0.0,0.0,(p.mpivot+p.mp)*p.g)
    accel_trans = M_trans \ B_trans
    x_ddot, y_ddot, z_ddot = accel_trans

    # -- Rotational dynamics about pivot --
    # Torque from bob forces about pivot
    lever = p.l * n
    tau = cross(lever, F .- Dp .- drag_bob)

    # Moment of inertia (point mass at distance l)
    I_theta = p.mp * p.l^2
    omega_dot = tau / I_theta
    omega_dot_x, omega_dot_y, omega_dot_z = omega_dot

    # -- Quaternion kinematics --
    # q_dot = ½ * Ω(omega) * q_vec
    Omega = @SMatrix [
         0.0      -omega_x   -omega_y   -omega_z;
      omega_x       0.0       omega_z   -omega_y;
      omega_y    -omega_z       0.0      omega_x;
      omega_z     omega_y     -omega_x      0.0
    ]
    q_vec = SVector(q0, q1, q2, q3)
    q_dot_vec = 0.5 * Omega * q_vec
    q0_dot, q1_dot, q2_dot, q3_dot = q_dot_vec

    # -- Pack derivatives --
    du[1]  = x_dot
    du[2]  = y_dot
    du[3]  = z_dot
    du[4]  = q0_dot
    du[5]  = q1_dot
    du[6]  = q2_dot
    du[7]  = q3_dot
    du[8]  = x_ddot
    du[9]  = y_ddot
    du[10] = z_ddot
    du[11] = omega_dot_x
    du[12] = omega_dot_y
    du[13] = omega_dot_z
end

# ----------------------- Set up & Solve ODE -------------------------
# Initial: upright quaternion [1,0,0,0], zero ang. vel.
u0 = [
    0.0,  0.0,  0.0,       # x,y,z
    1.0,  0.0,  0.0,  0.0, # q0,q1,q2,q3
    0.0,  0.0,  0.0,       # x_dot,y_dot,z_dot
    0.0,  0.0,  0.0        # omega_x,omega_y,omega_z
]
tspan = (0.0, 30.0)
prob = ODEProblem(pendulum_quat!, u0, tspan, params)

# Use a stiff method
sol = solve(prob, QNDF(), abstol=1e-6, reltol=1e-6)

println("Solution length: ", length(sol.u))

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
