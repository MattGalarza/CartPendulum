using GLMakie, Observables, DifferentialEquations, LinearAlgebra, Sundials, OrdinaryDiffEq, StaticArrays

# ----------------------------- System Parameters -----------------------------
struct Parameters
    mpivot::Float64        # pivot (rocket) mass
    mp::Float64            # pendulum mass
    l::Float64             # rod length
    g::Float64             # gravitational acceleration
    Bpivot::Float64        # pivot drag coefficient
    Bp::Float64            # pendulum drag coefficient
    control::Function      # (t,u) -> Vector{3}([Fx, Fy, Fz])
    disturbance::Function  # (t,u) -> Vector{3}([Dp_x, Dp_y, Dp_z])
end

zero_ctrl(t,u) = zeros(eltype(u),3)
zero_disturbance(t,u) = zeros(eltype(u),3)

params = Parameters(
    5.0,    # mpivot
    5.0,    # mp
    1.0,    # l
    9.81,   # g
    1.0,    # Bpivot
    1.0,    # Bp
    zero_ctrl,
    zero_disturbance
)

# ----------------------------- Quaternion Utilities -----------------------------
function euler_to_quaternion(θ, φ)
    θ2, φ2 = θ/2, φ/2
    qw =  cos(φ2)*cos(θ2)
    qx =  sin(φ2)*cos(θ2)
    qy =  sin(φ2)*sin(θ2)
    qz =  cos(φ2)*sin(θ2)
    return SVector(qw,qx,qy,qz)
end

function quaternion_to_direction(q::SVector{4,Float64})
    qw,qx,qy,qz = q
    dir_x = 2*(qx*qz + qw*qy)
    dir_y = 2*(qy*qz - qw*qx)
    dir_z = 1 - 2*(qx^2 + qy^2)
    return SVector(dir_x,dir_y,dir_z)
end

normalize_quaternion(q) = q / norm(q)

function quaternion_to_angles(q::SVector{4,Float64}, ω::AbstractVector)
    dir = quaternion_to_direction(q)
    φ = acos(clamp(dir[3], -1, 1))
    θ = atan(dir[2], dir[1]); θ < 0 && (θ += 2π)
    sφ = sin(φ)
    θ_dot = abs(sφ) < 1e-6 ? 0.0 : (ω[1]*cos(θ)+ω[2]*sin(θ))/sφ
    φ_dot = ω[1]*sin(θ) - ω[2]*cos(θ)
    return θ, φ, θ_dot, φ_dot
end

# ----------------------------- Dynamics Function -----------------------------
const ε = 1e-8  # regularization floor

function pendulum_quaternion!(du, u, p, t)
    x, y, z       = u[1:3]
    q             = normalize_quaternion(SVector(u[4:7]...))
    x_dot, y_dot, z_dot = u[8:10]
    ω             = SVector(u[11:13]...)

    # rod direction & lever arm
    dir   = quaternion_to_direction(q)
    r_vec = p.l .* dir

    # bob velocity = pivot vel + ω×r
    ω_cross_r = SVector(
      ω[2]*r_vec[3] - ω[3]*r_vec[2],
      ω[3]*r_vec[1] - ω[1]*r_vec[3],
      ω[1]*r_vec[2] - ω[2]*r_vec[1]
    )
    xp_dot, yp_dot, zp_dot = x_dot + ω_cross_r[1], y_dot + ω_cross_r[2], z_dot + ω_cross_r[3]

    # control & disturbance (adapter uses spherical for now)
    θ, φ, _, _ = quaternion_to_angles(q, ω)
    u_adapt    = [x,y,z,θ,φ,x_dot,y_dot,z_dot,0.0,0.0]
    F          = p.control(t, u_adapt)
    Dp         = p.disturbance(t, u_adapt)

    # drag
    drag_bob = p.Bp .* SVector(sign(xp_dot)*xp_dot^2,
                               sign(yp_dot)*yp_dot^2,
                               sign(zp_dot)*zp_dot^2)
    drag_piv = p.Bpivot .* SVector(sign(x_dot)*x_dot^2,
                                   sign(y_dot)*y_dot^2,
                                   sign(z_dot)*z_dot^2)

    gravity_force = SVector(0.0, 0.0, -p.mp * p.g)

    # --- rotational dynamics ---
    I3 = Matrix{Float64}(I,3,3)
    # build inertia_tensor in one shot (no in-place mutation)
    inertia_tensor = p.mp*p.l^2*(I3 - dir*dir') + ε*I3

    F_bob  = -gravity_force .- drag_bob .+ Dp
    torque = cross(r_vec, F_bob)
    L      = inertia_tensor * ω
    gyro   = cross(ω, L)
    angular_accel = inertia_tensor \ (torque - gyro)

    # --- translational coupling ---
    centripetal    = cross(ω, ω_cross_r)
    tangential     = cross(angular_accel, r_vec)
    reaction_accel = p.mp .* (centripetal .+ tangential)

    m_tot    = p.mpivot + p.mp
    pivot_acc = (SVector(F...) .- drag_piv .+ reaction_accel .- SVector(Dp...)) ./ m_tot
    pivot_acc   = pivot_acc - SVector(0.0,0.0,p.g)

    # --- quaternion kinematics ---
    Ω = @SMatrix [
       0.0   -ω[1]  -ω[2]  -ω[3];
      ω[1]    0.0    ω[3]  -ω[2];
      ω[2]   -ω[3]   0.0    ω[1];
      ω[3]    ω[2]  -ω[1]   0.0
    ]
    q_dot = 0.5 * (Ω * q)

    # --- pack derivatives elementwise ---
    du[1]  = x_dot
    du[2]  = y_dot
    du[3]  = z_dot
    du[4]  = q_dot[1]
    du[5]  = q_dot[2]
    du[6]  = q_dot[3]
    du[7]  = q_dot[4]
    du[8]  = pivot_acc[1]
    du[9]  = pivot_acc[2]
    du[10] = pivot_acc[3]
    du[11] = angular_accel[1]
    du[12] = angular_accel[2]
    du[13] = angular_accel[3]
end

# ----------------------------- Set up & Solve ODE Problem -----------------------------
# Convert initial spherical coordinates to quaternion
q0  = euler_to_quaternion(0.1, 0.3)

# Initial condition: [x, y, z, qw, qx, qy, qz, x_dot, y_dot, z_dot, ω_x, ω_y, ω_z]
z0 = vcat([0.0,0.0,0.0], q0, zeros(6))
tspan = (0.0, 35.0)

# Solve the pendulum system with increased data points
prob = ODEProblem(pendulum_quaternion!, z0, tspan, params)
sol  = solve(prob, CVODE_BDF(), abstol=1e-6, reltol=1e-6, maxiters=1000000, saveat=0.01)

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------- Create Figure + Plots -----------------------------
# Create 3D visualization
function animate_pendulum1(sol, params)
    # Create 3D visualization
    fig1 = Figure(size=(800, 800), fontsize=12)
    ax = fig1[1, 1] = Axis3(fig1, 
                         xlabel = "x", ylabel = "y", zlabel = "z",
                         limits = (-10*params.l, 10*params.l, -10*params.l, 10*params.l, -10*params.l, 10*params.l),
                         aspect = :data)
    
    # Get initial state
    u0 = sol.u[1]
    x, y, z = u0[1:3]
    
    # Convert to SVector explicitly
    q = SVector{4, Float64}(u0[4], u0[5], u0[6], u0[7])
    q = normalize_quaternion(q)
    
    # Get direction vector from quaternion
    dir = quaternion_to_direction(q)
    
    # Calculate pendulum position
    x_pend = x + params.l * dir[1]
    y_pend = y + params.l * dir[2]
    z_pend = z + params.l * dir[3]
    
    # Create visualization elements
    rocket_plot = meshscatter!(ax, [x], [y], [z], markersize = 0.2, color = :red)
    rod_plot = lines!(ax, [x, x_pend], [y, y_pend], [z, z_pend], linewidth = 3, color = :blue)
    pendulum_plot = meshscatter!(ax, [x_pend], [y_pend], [z_pend], markersize = 0.15, color = :blue)
    
    # Add quaternion visualization text
    quat_text = text!(ax, ["q = [" * join(round.([q[1], q[2], q[3], q[4]], digits=3), ", ") * "]"],
                    position = [(-2.5*params.l, -2.5*params.l, -2.5*params.l)],
                    color = :black, fontsize = 14)
    
    # Create trajectory visualization
    fig2 = Figure(size=(1200, 800))
    ax_3d_traj = fig2[1, 1] = Axis3(fig2, 
                                xlabel = "x", ylabel = "y", zlabel = "z",
                                title = "3D Trajectory",
                                limits = (-10*params.l, 10*params.l, -10*params.l, 10*params.l, -10*params.l, 10*params.l),
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
    
    # Initialize trajectory arrays within this function's scope
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
    
    # Display both figures
    display(GLMakie.Screen(), fig1)
    display(GLMakie.Screen(), fig2)
    
    # Animation parameters
    fps = 60
    dt_frame = 1/fps
    t_end = sol.t[end]
    
    # Animation loop
    sleep(1.0) # Brief delay to ensure windows are ready
    @async begin
        t_sim = 0.0
        
        while t_sim <= t_end && t_sim <= sol.t[end]
            try
                # Sample the solution at the current simulation timestep
                u = sol(t_sim)
                
                # Extract state components
                x, y, z = u[1:3]
                
                # Convert to SVector explicitly for quaternion
                q = SVector{4, Float64}(u[4], u[5], u[6], u[7])
                q = normalize_quaternion(q)
                
                x_dot, y_dot, z_dot = u[8:10]
                ω = SVector{3, Float64}(u[11], u[12], u[13])
                
                # Get direction vector from quaternion
                dir = quaternion_to_direction(q)
                
                # Calculate pendulum position
                x_pend = x + params.l * dir[1]
                y_pend = y + params.l * dir[2]
                z_pend = z + params.l * dir[3]
                
                # Convert quaternion to equivalent angles for phase plots
                theta_val, phi_val, theta_dot_val, phi_dot_val = quaternion_to_angles(q, ω)
                
                # Update visualization elements
                rocket_plot[1] = [x]
                rocket_plot[2] = [y]
                rocket_plot[3] = [z]
                
                rod_plot[1] = [x, x_pend]
                rod_plot[2] = [y, y_pend]
                rod_plot[3] = [z, z_pend]
                
                pendulum_plot[1] = [x_pend]
                pendulum_plot[2] = [y_pend]
                pendulum_plot[3] = [z_pend]
                
                # Update quaternion text display
                quat_text[1] = ["q = [" * join(round.([q[1], q[2], q[3], q[4]], digits=3), ", ") * "]"]
                
                # Store trajectory data
                push!(rocket_x, x)
                push!(rocket_y, y)
                push!(rocket_z, z)
                push!(pendulum_x, x_pend)
                push!(pendulum_y, y_pend)
                push!(pendulum_z, z_pend)
                
                # Store phase portrait data
                push!(theta, theta_val)
                push!(thetadot, theta_dot_val)
                push!(phi, phi_val)
                push!(phidot, phi_dot_val)
                push!(time_array, t_sim)
                push!(z_height, z)
                
                # Keep trajectories a reasonable length and ensure all arrays have same length
                max_points = length(sol.u)
                
                if length(rocket_x) > max_points
                    # Calculate minimum length to ensure all arrays have same size
                    min_length = min(
                        length(rocket_x), length(rocket_y), length(rocket_z),
                        length(pendulum_x), length(pendulum_y), length(pendulum_z),
                        length(theta), length(thetadot),
                        length(phi), length(phidot),
                        length(time_array), length(z_height)
                    )
                    
                    # Ensure min_length doesn't exceed max_points
                    min_length = min(min_length, max_points)
                    
                    # Trim all arrays to the same length
                    rocket_x = rocket_x[end-min_length+1:end]
                    rocket_y = rocket_y[end-min_length+1:end]
                    rocket_z = rocket_z[end-min_length+1:end]
                    pendulum_x = pendulum_x[end-min_length+1:end]
                    pendulum_y = pendulum_y[end-min_length+1:end]
                    pendulum_z = pendulum_z[end-min_length+1:end]
                    theta = theta[end-min_length+1:end]
                    thetadot = thetadot[end-min_length+1:end]
                    phi = phi[end-min_length+1:end]
                    phidot = phidot[end-min_length+1:end]
                    time_array = time_array[end-min_length+1:end]
                    z_height = z_height[end-min_length+1:end]
                end

                # Update trajectory and phase plots
                rocket_traj[1] = rocket_x
                rocket_traj[2] = rocket_y
                rocket_traj[3] = rocket_z
                pendulum_traj[1] = pendulum_x
                pendulum_traj[2] = pendulum_y
                pendulum_traj[3] = pendulum_z
                
                theta_line[1] = theta
                theta_line[2] = thetadot
                phi_line[1] = phi
                phi_line[2] = phidot
                z_line[1] = time_array
                z_line[2] = z_height

                sleep(dt_frame)
                t_sim += dt_frame
            catch e
                println("Error at t=$t_sim: $e")
                println("Error type: ", typeof(e))
                println("Full error: ", sprint(showerror, e))
                break
            end
        end
    end
    
    println("3D Pendulum simulation is running!")
    return fig1, fig2  # Return the figures so they stay alive
end

function animate_pendulum2(sol, params)
    # Create 3D visualization
    fig1 = Figure(size=(800, 800), fontsize=12)
    ax = fig1[1, 1] = Axis3(fig1, 
                         xlabel = "x", ylabel = "y", zlabel = "z",
                         limits = (-10*params.l, 10*params.l, -10*params.l, 10*params.l, -10*params.l, 10*params.l),
                         aspect = :data)

    # Get initial state
    u0 = sol.u[1]
    x, y, z = u0[1:3]

    q = normalize_quaternion(SVector{4, Float64}(u0[4:7]))
    dir = quaternion_to_direction(q)

    # Calculate pendulum position
    x_pend = x + params.l * dir[1]
    y_pend = y + params.l * dir[2]
    z_pend = z + params.l * dir[3]

    # Create visualization elements
    rocket_plot = meshscatter!(ax, [x], [y], [z], markersize = 0.2, color = :red)
    rod_plot = lines!(ax, [x, x_pend], [y, y_pend], [z, z_pend], linewidth = 3, color = :blue)
    pendulum_plot = meshscatter!(ax, [x_pend], [y_pend], [z_pend], markersize = 0.15, color = :blue)

    quat_text = text!(ax, ["q = [" * join(round.([q[1], q[2], q[3], q[4]], digits=3), ", ") * "]"],
                    position = [(-2.5*params.l, -2.5*params.l, -2.5*params.l)],
                    color = :black, fontsize = 14)

    # Create trajectory visualization
    fig2 = Figure(size=(1200, 800))
    ax_3d_traj = fig2[1:2, 1] = Axis3(fig2, 
                                xlabel = "x", ylabel = "y", zlabel = "z",
                                title = "3D Trajectory",
                                limits = (-10*params.l, 10*params.l, -10*params.l, 10*params.l, -10*params.l, 10*params.l),
                                aspect = :data)

    ax_empty1 = fig2[1, 2] = Axis(fig2, title="Empty Plot 1")
    ax_empty2 = fig2[2, 2] = Axis(fig2, title="Empty Plot 2")

    rocket_traj = lines!(ax_3d_traj, Float64[], Float64[], Float64[], color = :red, label = "Pivot")
    pendulum_traj = lines!(ax_3d_traj, Float64[], Float64[], Float64[], color = :blue, label = "Pendulum")
    Legend(fig2[1, 3], ax_3d_traj)

    rocket_x, rocket_y, rocket_z = Float64[], Float64[], Float64[]
    pendulum_x, pendulum_y, pendulum_z = Float64[], Float64[], Float64[]

    display(GLMakie.Screen(), fig1)
    display(GLMakie.Screen(), fig2)

    fps, dt_frame, t_end = 60, 1/60, sol.t[end]

    sleep(1.0)
    @async begin
        t_sim = 0.0
        while t_sim <= t_end && t_sim <= sol.t[end]
            try
                u = sol(t_sim)

                x, y, z = u[1:3]
                q = normalize_quaternion(SVector{4, Float64}(u[4:7]))
                dir = quaternion_to_direction(q)

                x_pend = x + params.l * dir[1]
                y_pend = y + params.l * dir[2]
                z_pend = z + params.l * dir[3]

                rocket_plot[1] = [x]
                rocket_plot[2] = [y]
                rocket_plot[3] = [z]

                rod_plot[1] = [x, x_pend]
                rod_plot[2] = [y, y_pend]
                rod_plot[3] = [z, z_pend]

                pendulum_plot[1] = [x_pend]
                pendulum_plot[2] = [y_pend]
                pendulum_plot[3] = [z_pend]

                quat_text[1] = ["q = [" * join(round.([q[1], q[2], q[3], q[4]], digits=3), ", ") * "]"]

                push!(rocket_x, x)
                push!(rocket_y, y)
                push!(rocket_z, z)
                push!(pendulum_x, x_pend)
                push!(pendulum_y, y_pend)
                push!(pendulum_z, z_pend)

                rocket_traj[1] = rocket_x
                rocket_traj[2] = rocket_y
                rocket_traj[3] = rocket_z

                pendulum_traj[1] = pendulum_x
                pendulum_traj[2] = pendulum_y
                pendulum_traj[3] = pendulum_z

                ax_3d_traj.limits = (x-5params.l, x+5params.l, y-5params.l, y+5params.l, z-5params.l, z+5params.l)
                ax.limits = (x-5params.l, x+5params.l, y-5params.l, y+5params.l, z-5params.l, z+5params.l)

                sleep(dt_frame)
                t_sim += dt_frame
            catch e
                println("Error at t=$t_sim: $e")
                break
            end
        end
    end

    println("3D Pendulum simulation is running!")
    return fig1, fig2
end

# After solving the ODE, call the animation function
println("Setting up animation...")
figs = animate_pendulum2(sol, params)




















# ----------------------------- System Parameters -----------------------------
struct Parameters
    mpivot::Float64        # pivot (rocket) mass
    mp::Float64            # pendulum mass
    l::Float64             # rod length
    g::Float64             # gravitational acceleration
    Bpivot::Float64        # pivot drag coefficient
    Bp::Float64            # pendulum drag coefficient
    control::Function      # (t,u) -> Vector{3}([Fx, Fy, Fz])
    disturbance::Function  # (t,u) -> Vector{3}([Dp_x, Dp_y, Dp_z])
end

zero_ctrl(t,u) = zeros(eltype(u),3)
zero_disturbance(t,u) = zeros(eltype(u),3)

params = Parameters(
    5.0,    # mpivot
    5.0,    # mp
    1.0,    # l
    9.81,   # g
    1.0,    # Bpivot
    1.0,    # Bp
    zero_ctrl,
    zero_disturbance
)

# ----------------------------- Quaternion Utilities -----------------------------
function euler_to_quaternion(θ, φ)
    θ2, φ2 = θ/2, φ/2
    qw =  cos(φ2)*cos(θ2)
    qx =  sin(φ2)*cos(θ2)
    qy =  sin(φ2)*sin(θ2)
    qz =  cos(φ2)*sin(θ2)
    return SVector(qw,qx,qy,qz)
end

function quaternion_to_direction(q::SVector{4,Float64})
    qw,qx,qy,qz = q
    dir_x = 2*(qx*qz + qw*qy)
    dir_y = 2*(qy*qz - qw*qx)
    dir_z = 1 - 2*(qx^2 + qy^2)
    return SVector(dir_x,dir_y,dir_z)
end

normalize_quaternion(q) = q / norm(q)

function quaternion_to_angles(q::SVector{4,Float64}, ω::SVector{3,Float64})
    dir = quaternion_to_direction(q)
    φ = acos(clamp(dir[3], -1, 1))
    θ = atan(dir[2], dir[1]); θ < 0 && (θ += 2π)
    sφ = sin(φ)
    θ_dot = abs(sφ) < 1e-6 ? 0.0 : (ω[1]*cos(θ) + ω[2]*sin(θ)) / sφ
    φ_dot = ω[1]*sin(θ) - ω[2]*cos(θ)
    return θ, φ, θ_dot, φ_dot
end

# ----------------------------- Dynamics Function -----------------------------
const ε = 1e-8  # regularization floor

function pendulum_quaternion!(du, u, p, t)
    # Unpack states
    x, y, z       = u[1:3]
    q             = normalize_quaternion(SVector(u[4:7]...))
    x_dot, y_dot, z_dot = u[8:10]
    ω             = SVector(u[11:13]...)

    # rod direction & lever arm
    dir   = quaternion_to_direction(q)
    r_vec = p.l .* dir

    # bob velocity = pivot vel + ω × r
    ω_cross_r = SVector(
      ω[2]*r_vec[3] - ω[3]*r_vec[2],
      ω[3]*r_vec[1] - ω[1]*r_vec[3],
      ω[1]*r_vec[2] - ω[2]*r_vec[1]
    )
    xp_dot = x_dot + ω_cross_r[1]
    yp_dot = y_dot + ω_cross_r[2]
    zp_dot = z_dot + ω_cross_r[3]

    # spherical‑adapter now includes θ̇, φ̇
    θ, φ, θ_dot, φ_dot = quaternion_to_angles(q, ω)
    u_adapt    = [x, y, z, θ, φ, x_dot, y_dot, z_dot, θ_dot, φ_dot]

    # control & disturbance
    F   = p.control(t, u_adapt)
    Dp  = p.disturbance(t, u_adapt)

    # drag forces
    drag_bob = p.Bp .* SVector(sign(xp_dot)*xp_dot^2,
                               sign(yp_dot)*yp_dot^2,
                               sign(zp_dot)*zp_dot^2)
    drag_piv = p.Bpivot .* SVector(sign(x_dot)*x_dot^2,
                                   sign(y_dot)*y_dot^2,
                                   sign(z_dot)*z_dot^2)

    # gravity on bob
    gravity_force = SVector(0.0, 0.0, -p.mp * p.g)

    # net force on bob → torque
    F_bob  = gravity_force .- drag_bob .+ Dp

    # rotational dynamics
    I3s = @SMatrix [1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 0.0 1.0]
    inertia_tensor = p.mp * p.l^2 * (I3s - dir*dir') + ε * I3s

    torque  = cross(r_vec, F_bob)
    L       = inertia_tensor * ω
    gyro    = cross(ω, L)
    α       = inertia_tensor \ (torque - gyro)  # angular accel

    # translational coupling
    centripetal    = cross(ω, ω_cross_r)
    tangential     = cross(α, r_vec)
    reaction_accel = p.mp .* (centripetal .+ tangential)

    # pivot acceleration (no double‑subtraction of Dp)
    m_tot    = p.mpivot + p.mp
    pivot_acc = (SVector(F...) .- drag_piv .+ reaction_accel) ./ m_tot
    pivot_acc -= SVector(0.0, 0.0, p.g)

    # quaternion kinematics
    Ω = @SMatrix [
       0.0   -ω[1]  -ω[2]  -ω[3];
      ω[1]    0.0    ω[3]  -ω[2];
      ω[2]   -ω[3]   0.0    ω[1];
      ω[3]    ω[2]  -ω[1]   0.0
    ]
    q_dot = 0.5 * (Ω * q)

    # pack derivatives
    du[1:3]   .= (x_dot, y_dot, z_dot)
    du[4:7]   .= q_dot
    du[8:10]  .= pivot_acc
    du[11:13] .= α
end

# ----------------------------- Quaternion Renormalization Callback -------------
function renormalize_q!(integrator)
    u = integrator.u
    q = SVector{4,Float64}(u[4], u[5], u[6], u[7])
    qn = q / norm(q)
    integrator.u[4:7] = qn
end
cb = DiscreteCallback((u,t,integrator) -> true, renormalize_q!)

# ----------------------------- Set up & Solve ODE Problem -----------------------------
# initial quaternion from small Euler tilt
q0  = euler_to_quaternion(0.1, 0.3)

# state: [x,y,z, qw,qx,qy,qz, ẋ,ẏ,ż, ωx,ωy,ωz]
z0 = vcat([0.0, 0.0, 0.0], q0, zeros(6))
tspan = (0.0, 35.0)

prob = ODEProblem(pendulum_quaternion!, z0, tspan, params)

sol = solve(
    prob,
    CVODE_BDF(),
    callback = cb,
    abstol   = 1e-6,
    reltol   = 1e-6,
    maxiters = 1000000,
    saveat   = 0.01
)

println("Solver returned: ", sol.retcode)
