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
    0.3,    # Bpivot
    2.0,    # Bp
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
    x, y, z = u[1:3]
    q = normalize_quaternion(SVector(u[4:7]...))
    x_dot, y_dot, z_dot = u[8:10]
    ω = SVector(u[11:13]...)

    # rod direction & lever arm
    dir = quaternion_to_direction(q)
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
    u_adapt = [x,y,z,θ,φ,x_dot,y_dot,z_dot,0.0,0.0]
    F = p.control(t, u_adapt)
    Dp = p.disturbance(t, u_adapt)

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

    F_bob = gravity_force .- drag_bob .+ Dp
    torque = cross(r_vec, F_bob)
    L = inertia_tensor * ω
    gyro = cross(ω, L)
    angular_accel = inertia_tensor \ (torque - gyro)

    # --- translational coupling ---
    centripetal = cross(ω, ω_cross_r)
    tangential = cross(angular_accel, r_vec)
    reaction_accel = p.mp .* (centripetal .+ tangential)

    m_tot = p.mpivot + p.mp
    pivot_acc = (SVector(F...) .- drag_piv .+ reaction_accel .- SVector(Dp...)) ./ m_tot
    pivot_acc = pivot_acc - SVector(0.0,0.0,p.g)

    # --- quaternion kinematics ---
    Ω = @SMatrix [
       0.0   -ω[1]  -ω[2]  -ω[3];
      ω[1]    0.0    ω[3]  -ω[2];
      ω[2]   -ω[3]   0.0    ω[1];
      ω[3]    ω[2]  -ω[1]   0.0
    ]
    q_dot = 0.5 * (Ω * q)

    # --- pack derivatives elementwise ---
    du[1] = x_dot
    du[2] = y_dot
    du[3] = z_dot
    du[4] = q_dot[1]
    du[5] = q_dot[2]
    du[6] = q_dot[3]
    du[7] = q_dot[4]
    du[8] = pivot_acc[1]
    du[9] = pivot_acc[2]
    du[10] = pivot_acc[3]
    du[11] = angular_accel[1]
    du[12] = angular_accel[2]
    du[13] = angular_accel[3]
end

# ----------------------------- Set up & Solve ODE Problem -----------------------------
# Convert initial spherical coordinates to quaternion
q0 = euler_to_quaternion(0.1, 0.3)

# Initial condition: [x, y, z, qw, qx, qy, qz, x_dot, y_dot, z_dot, ω_x, ω_y, ω_z]
z0 = vcat([0.0,0.0,0.0], q0, zeros(6))
tspan = (0.0, 35.0)

# Solve the pendulum system with increased data points
prob = ODEProblem(pendulum_quaternion!, z0, tspan, params)
sol = solve(prob, CVODE_BDF(), abstol=1e-6, reltol=1e-6, maxiters=1000000, saveat=0.01)

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

# pick the figure you animated:
fig = fig1  

# define the times you want to capture (0…t_end in N frames)
t_end = sol.t[end]
nframes = Int(round(60 * t_end))
times = range(0, stop=t_end, length=nframes)

# record to GIF
record(fig, "pendulum.gif", times) do t
    sol_t = sol(t)
    # update your scene exactly as you do inside your @async loop:
    # e.g. compute x,y,z,q,…, update meshscatter!, lines!, text!, etc.
end

# or record to MP4 (requires ffmpeg on your PATH):
record(fig, "pendulum.mp4", times; framerate=60)





using GLMakie, Observables, DifferentialEquations, LinearAlgebra, Sundials, OrdinaryDiffEq, StaticArrays

# ----------------------------- System Parameters -----------------------------
struct Parameters
    mpivot::Float64        # pivot (rocket) mass
    mp::Float64            # pendulum mass
    l::Float64             # rod length
    g::Float64             # gravitational acceleration
    Bpivot::Float64        # pivot drag coefficient
    Bp::Float64            # pendulum drag coefficient
    Brot::Float64          # rotational damping
    Kp::Float64            # proportional gain
    Ki::Float64            # integral gain
    Kd::Float64            # derivative gain
    disturbance::Function  # (t,u) -> Vector{3}([Dp_x, Dp_y, Dp_z])
end

zero_disturbance(t,u) = SVector{3,Float64}(0.0,0.0,0.0)

params = Parameters(
    5.0,    # mpivot
    5.0,    # mp
    1.0,    # l
    9.81,   # g
    5.0,    # Bpivot
    5.0,    # Bp
    0.3,    # Brot
    0.0,    # Kp
    0.0,    # Ki
    0.0,   # Kd
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

# ----------------------------- Dynamics + PID Control -----------------------------
const ε = 1e-8  # regularization floor

function pendulum_quaternion!(du, u, p, t)
    # unpack state (now 16 components)
    x, y, z = u[1:3]
    q = normalize_quaternion(SVector(u[4:7]...))
    x_dot, y_dot, z_dot = u[8:10]
    ω = SVector(u[11:13]...)
    σx, σy, σz = u[14:16]    # integrators

    # rod direction & lever arm
    dir = quaternion_to_direction(q)
    r_vec = p.l .* dir

    # bob velocity = pivot vel + ω×r
    ωxr = SVector(
      ω[2]*r_vec[3] - ω[3]*r_vec[2],
      ω[3]*r_vec[1] - ω[1]*r_vec[3],
      ω[1]*r_vec[2] - ω[2]*r_vec[1]
    )
    xp_dot = x_dot + ωxr[1]
    yp_dot = y_dot + ωxr[2]
    zp_dot = z_dot + ωxr[3]

    # position of bob
    x_p = x + p.l*dir[1]
    y_p = y + p.l*dir[2]
    z_p = z + p.l*dir[3]

    # Setpoint positions
    X = 0
    Y = 0
    Z = p.l

    # PID error
    ex = x_p - (x + X)  
    ey = y_p - (y + Y)  
    ez = z_p - (z + Z) 

    dex, dey, dez = xp_dot - x_dot, yp_dot - y_dot, zp_dot - z_dot

    Fx = -p.Kp*ex - p.Ki*σx - p.Kd*dex
    Fy = -p.Kp*ey - p.Ki*σy - p.Kd*dey
    Fz = -p.Kp*ez - p.Ki*σz - p.Kd*dez # (p.mp + p.mpivot)*p.g

    F = SVector(Fx, Fy, Fz)
    Dp = SVector(p.disturbance(t,u)...)

    # drag
    drag_bob = p.Bp .* SVector(sign(xp_dot)*xp_dot^2,
                               sign(yp_dot)*yp_dot^2,
                               sign(zp_dot)*zp_dot^2)
    drag_piv = p.Bpivot .* SVector(sign(x_dot)*x_dot^2,
                                   sign(y_dot)*y_dot^2,
                                   sign(z_dot)*z_dot^2)

    # --- rotational dynamics ---
    I3 = I*1.0 # identity 3×3
    inertia_tensor = p.mp*p.l^2*(I3 - dir*dir') + ε*I3

    # gravity on bob
    gravity_force = SVector(0.0, 0.0, -p.mp * p.g)

    # net force on bob
    F_bob = gravity_force .- drag_bob .+ Dp

    # rotational dynamics 
    torque = cross(r_vec, F_bob)
    L = inertia_tensor * ω
    gyro = cross(ω, L)
    torque_damped = torque .- p.Brot .* ω
    angular_accel = inertia_tensor \ (torque_damped .- gyro)

    # translational coupling
    centripetal = cross(ω, ωxr)
    tangential = cross(angular_accel, r_vec)

    # rod‐reaction on pivot is *negative* of bob’s internal acceleration
    reaction_force = -p.mp .* (centripetal .+ tangential)

    # pivot acceleration
    m_tot = p.mpivot + p.mp
    pivot_acc = (F .- drag_piv .+ reaction_force) ./ m_tot
    pivot_acc -= SVector(0.0, 0.0, p.g)

    # --- quaternion kinematics ---
    Ω    = @SMatrix [
       0.0   -ω[1]  -ω[2]  -ω[3];
      ω[1]    0.0    ω[3]  -ω[2];
      ω[2]   -ω[3]   0.0    ω[1];
      ω[3]    ω[2]  -ω[1]   0.0
    ]
    q_dot = 0.5 * (Ω * q)

    # --- pack derivatives ---
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
    du[14] = ex     # integral updates
    du[15] = ey
    du[16] = ez
end

# ----------------------------- Quaternion Renormalization Callback -------------
function renormalize_q!(integrator)
  u = integrator.u
  q = SVector{4,Float64}(u[4],u[5],u[6],u[7])
  integrator.u[4:7] = q/norm(q)
end
cb = DiscreteCallback((u,t,integrator)->true, renormalize_q!)

# ----------------------------- Set up & Solve ODE Problem -----------------------------
# Convert initial spherical coordinates to quaternion
q0 = euler_to_quaternion(0.0, 0.1)

# Initial condition
z0 = vcat([0.0,0.0,0.0], q0, zeros(6), [0.0,0.0,0.0])  
tspan = (0.0, 25.0)

prob = ODEProblem(pendulum_quaternion!, z0, tspan, params)
sol  = solve(prob, CVODE_BDF(), abstol=1e-6, reltol=1e-6, maxiters=10^6, saveat=0.01)

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

function animate_pendulum_combined(sol, params; save_video=true, filename="pendulum_combined2.mp4")
    # Create combined figure
    fig_combined = Figure(size=(1200, 600), fontsize=12)
    
    # Left panel: 3D pendulum visualization
    ax_pendulum = fig_combined[1, 1] = Axis3(fig_combined, 
                        xlabel = "x", ylabel = "y", zlabel = "z",
                        title = "Pendulum Motion",
                        limits = (-10*params.l, 10*params.l, -10*params.l, 10*params.l, -10*params.l, 10*params.l),
                        aspect = :data)
    
    # Right panel: 3D trajectory with FIXED window size but moving center
    ax_trajectory = fig_combined[1, 2] = Axis3(fig_combined, 
                        xlabel = "x", ylabel = "y", zlabel = "z",
                        title = "3D Trajectory",
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
    
    # Create visualization elements for pendulum
    rocket_plot = meshscatter!(ax_pendulum, [x], [y], [z], markersize = 0.2, color = :red)
    rod_plot = lines!(ax_pendulum, [x, x_pend], [y, y_pend], [z, z_pend], linewidth = 3, color = :blue)
    pendulum_plot = meshscatter!(ax_pendulum, [x_pend], [y_pend], [z_pend], markersize = 0.15, color = :blue)
    
    quat_text = text!(ax_pendulum, ["q = [" * join(round.([q[1], q[2], q[3], q[4]], digits=3), ", ") * "]"],
                position = [(-2.5*params.l, -2.5*params.l, -2.5*params.l)],
                color = :black, fontsize = 14)
    
    # Create trajectories for the plot
    rocket_traj = lines!(ax_trajectory, Float64[], Float64[], Float64[], color = :red, label = "Pivot")
    pendulum_traj = lines!(ax_trajectory, Float64[], Float64[], Float64[], color = :blue, label = "Pendulum")
    Legend(fig_combined[1, 3], ax_trajectory)
    
    rocket_x, rocket_y, rocket_z = Float64[], Float64[], Float64[]
    pendulum_x, pendulum_y, pendulum_z = Float64[], Float64[], Float64[]
    
    # Set up video recording
    fps, dt_frame, t_end = 60, 1/60, sol.t[end]
    frames = Observable(1)
    
    # Define the size of the viewing window for the trajectory plot
    window_size = 10 * params.l
    
    record(fig_combined, filename, 1:Int(ceil(t_end/dt_frame)); framerate = fps) do frame
        t_sim = (frame - 1) * dt_frame
        if t_sim <= t_end
            u = sol(t_sim)
            
            x, y, z = u[1:3]
            q = normalize_quaternion(SVector{4, Float64}(u[4:7]))
            dir = quaternion_to_direction(q)
            
            x_pend = x + params.l * dir[1]
            y_pend = y + params.l * dir[2]
            z_pend = z + params.l * dir[3]
            
            # Update pendulum visualization
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
            
            # Update trajectory
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
            
            # Update view limits for BOTH plots to follow their respective centers
            # For pendulum plot, keep it centered on current pendulum position
            ax_pendulum.limits = (x-5*params.l, x+5*params.l, y-5*params.l, y+5*params.l, z-5*params.l, z+5*params.l)
            
            # For trajectory plot, keep a fixed-size window centered on the current position
            ax_trajectory.limits = (x-window_size, x+window_size, 
                                    y-window_size, y+window_size, 
                                    z-window_size, z+window_size)
        end
        frames[] = frame
    end
    
    println("Combined video saved to $filename")
    return fig_combined
end

# After solving the ODE, call the animation function
println("Creating combined pendulum and trajectory video...")
fig_combined = animate_pendulum_combined(sol, params, filename="pendulum_combined2.mp4")

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




# 1) Build arrays of t, Fx, Fy, Fz
ts   = sol.t
Fx_a = Float64[]
Fy_a = Float64[]
Fz_a = Float64[]

for u in sol.u
    # reconstruct geometry exactly as in the integrator
    dir       = quaternion_to_direction(normalize_quaternion(SVector(u[4:7]...)))
    r_vec     = params.l * dir
    ω         = SVector(u[11], u[12], u[13])
    ωxr       = SVector(
                  ω[2]*r_vec[3] - ω[3]*r_vec[2],
                  ω[3]*r_vec[1] - ω[1]*r_vec[3],
                  ω[1]*r_vec[2] - ω[2]*r_vec[1]
                )
    # positions & velocities of bob
    xp        = u[1] + r_vec[1]
    yp        = u[2] + r_vec[2]
    zp        = u[3] + r_vec[3]
    xp_dot    = u[8]  + ωxr[1]
    yp_dot    = u[9]  + ωxr[2]
    zp_dot    = u[10] + ωxr[3]

    # errors & derivatives
    e_x,  e_y,  e_z  = xp - u[1], yp - u[2], zp - u[3]
    de_x, de_y, de_z = xp_dot - u[8], yp_dot - u[9], zp_dot - u[10]

    # PID formula
    Fx = -params.Kp*e_x - params.Ki*u[14] - params.Kd*de_x
    Fy = -params.Kp*e_y - params.Ki*u[15] - params.Kd*de_y
    Fz = -params.Kp*e_z - params.Ki*u[16] - params.Kd*de_z

    push!(Fx_a, Fx)
    push!(Fy_a, Fy)
    push!(Fz_a, Fz)
end

fig = Figure(resolution = (800,600))
ax  = Axis(fig[1,1],
           xlabel = "Time (s)",
           ylabel = "Force (N)",
           title  = "PID Control Inputs Over Time")

lines!(ax, ts, Fx_a, label = "Fₓ")
lines!(ax, ts, Fy_a, label = "Fᵧ")
lines!(ax, ts, Fz_a, label = "F_z")

axislegend(ax)
display(fig)

function plot_pendulum_static(sol, p)
    # reference height for z‑error (initial pivot z + rod length)
    z0    = sol.u[1][3]
    z_ref = z0 + p.l

    ts  = sol.t
    N   = length(ts)

    # preallocate arrays
    px  = Vector{Float64}(undef, N)
    py  = Vector{Float64}(undef, N)
    pz  = Vector{Float64}(undef, N)
    bx  = Vector{Float64}(undef, N)
    by  = Vector{Float64}(undef, N)
    bz  = Vector{Float64}(undef, N)
    Fx  = Vector{Float64}(undef, N)
    Fy  = Vector{Float64}(undef, N)
    Fz  = Vector{Float64}(undef, N)
    ex  = Vector{Float64}(undef, N)
    ey  = Vector{Float64}(undef, N)
    ez  = Vector{Float64}(undef, N)

    for (i, t) in enumerate(ts)
        u = sol(t)

        # unpack pivot
        x, y, z = u[1], u[2], u[3]

        # quaternion → direction
        q   = normalize_quaternion(SVector(u[4],u[5],u[6],u[7]))
        dir = quaternion_to_direction(q)

        # bob position
        x_p = x + p.l * dir[1]
        y_p = y + p.l * dir[2]
        z_p = z + p.l * dir[3]

        # bob velocity = pivot vel + ω×r
        ω     = SVector(u[11],u[12],u[13])
        r_vec = p.l * dir
        ωxr   = SVector(
            ω[2]*r_vec[3] - ω[3]*r_vec[2],
            ω[3]*r_vec[1] - ω[1]*r_vec[3],
            ω[1]*r_vec[2] - ω[2]*r_vec[1]
        )
        xp_dot = u[8]  + ωxr[1]
        yp_dot = u[9]  + ωxr[2]
        zp_dot = u[10] + ωxr[3]

        # errors for PID
        ex[i] = x_p - x
        ey[i] = y_p - y
        ez[i] = z_ref - z_p

        # derivative of error
        dex = xp_dot - u[8]
        dey = yp_dot - u[9]
        dez = zp_dot - u[10]

        # integrator states
        σx, σy, σz = u[14], u[15], u[16]

        # PID outputs
        Fx[i] = -p.Kp*ex[i] - p.Ki*σx - p.Kd*dex
        Fy[i] = -p.Kp*ey[i] - p.Ki*σy - p.Kd*dey
        Fz[i] = -p.Kp*ez[i] - p.Ki*σz - p.Kd*dez

        # store trajectories
        px[i], py[i], pz[i] = x, y, z
        bx[i], by[i], bz[i] = x_p, y_p, z_p
    end

    # now build the figure
    fig = Figure(resolution=(1200,800))

    # 3D trajectories
    ax3 = Axis3(fig[1:2,1], 
        xlabel="x", ylabel="y", zlabel="z", title="3D Trajectories", aspect=:data)
    lines!(ax3, px, py, pz, color=:red,  label="Pivot")
    lines!(ax3, bx, by, bz, color=:blue, label="Bob")
    Legend(fig[1,2], ax3)

    # PID forces vs time
    axF = Axis(fig[1,2],
        xlabel="t (s)", ylabel="F (N)", title="PID Control Forces")
    lines!(axF, ts, Fx, color=:blue,  label="Fₓ")
    lines!(axF, ts, Fy, color=:orange,label="Fᵧ")
    lines!(axF, ts, Fz, color=:green, label="F_z")
    axislegend(axF)

    # Position errors vs time
    axE = Axis(fig[2,2],
        xlabel="t (s)", ylabel="Error (m)", title="Position Errors")
    lines!(axE, ts, ex, color=:blue,  label="eₓ")
    lines!(axE, ts, ey, color=:orange,label="eᵧ")
    lines!(axE, ts, ez, color=:green, label="e_z")
    axislegend(axE)

    display(fig)
    return fig
end

# After solving the ODE:
plot_pendulum_static(sol, params)

function plot_errors(sol, p)
    # reference height for z‑error (initial pivot z + rod length)
    z0    = sol.u[1][3]
    z_ref = z0 + p.l

    ts  = sol.t
    N   = length(ts)

    # preallocate error arrays
    ex = Vector{Float64}(undef, N)
    ey = Vector{Float64}(undef, N)
    ez = Vector{Float64}(undef, N)

    for (i, t) in enumerate(ts)
        u = sol(t)

        # pivot pos
        x, y, z = u[1], u[2], u[3]

        # rod orientation → bob pos
        q   = normalize_quaternion(SVector(u[4],u[5],u[6],u[7]))
        dir = quaternion_to_direction(q)
        x_p = x + p.l*dir[1]
        y_p = y + p.l*dir[2]
        z_p = z + p.l*dir[3]

        # errors
        ex[i] = x_p - x
        ey[i] = y_p - y
        ez[i] = z_ref - z_p
    end

    fig = Figure(resolution=(800,400))
    ax  = Axis(fig[1,1], xlabel="t (s)", ylabel="Error (m)", title="Position Errors")
    lines!(ax, ts, ex, color=:blue,  label="eₓ")
    lines!(ax, ts, ey, color=:orange,label="eᵧ")
    lines!(ax, ts, ez, color=:green, label="e_z")
    axislegend(ax)

    display(fig)
    return fig
end

# After solving:
plot_errors(sol, params)







using GLMakie, Observables, DifferentialEquations, LinearAlgebra, Sundials, OrdinaryDiffEq, StaticArrays

# ----------------------------- System Parameters -----------------------------
struct Parameters
    mpivot::Float64        # pivot (rocket) mass
    mp::Float64            # pendulum mass
    l::Float64             # rod length
    g::Float64             # gravitational acceleration
    Bpivot::Float64        # pivot drag coefficient
    Bp::Float64            # pendulum drag coefficient
    Brot::Float64          # rotational damping
    Kp::Float64            # proportional gain
    Ki::Float64            # integral gain
    Kd::Float64            # derivative gain
    disturbance::Function  # (t,u) -> Vector{3}([Dp_x, Dp_y, Dp_z])
end

zero_disturbance(t,u) = SVector{3,Float64}(0.0, 0.0, 0.0)

params = Parameters(
    5.0,    # mpivot
    5.0,    # mp
    1.0,    # l
    9.81,   # g
    5.0,    # Bpivot
    5.0,    # Bp
    1.5,    # Brot
    120.0,  # Kp
    1.0,    # Ki
    40.0,   # Kd
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
    qw, qx, qy, qz = q
    dir_x = 2*(qx*qz + qw*qy)
    dir_y = 2*(qy*qz - qw*qx)
    dir_z = 1 - 2*(qx^2 + qy^2)
    return SVector(dir_x,dir_y,dir_z)
end

normalize_quaternion(q) = q / max(norm(q), 1e-10)  # Avoid division by near-zero

function quaternion_to_angles(q::SVector{4,Float64}, ω::AbstractVector)
    dir = quaternion_to_direction(q)
    φ = acos(clamp(dir[3], -1, 1))
    θ = atan(dir[2], dir[1]); θ < 0 && (θ += 2π)
    sφ = sin(φ)
    θ_dot = abs(sφ) < 1e-6 ? 0.0 : (ω[1]*cos(θ)+ω[2]*sin(θ))/sφ
    φ_dot = ω[1]*sin(θ) - ω[2]*cos(θ)
    return θ, φ, θ_dot, φ_dot
end

# ----------------------------- Dynamics + Energy-Based Control ------------------
const ε = 1e-10  # regularization floor

function pendulum_quaternion!(du, u, p, t)
    # unpack state (now 16 components)
    x, y, z = u[1:3]
    q = SVector{4,Float64}(u[4], u[5], u[6], u[7])
    x_dot, y_dot, z_dot = u[8:10]
    ω = SVector{3,Float64}(u[11], u[12], u[13])
    σx, σy, σz = u[14:16]    # integrators

    # rod direction & lever arm
    dir = quaternion_to_direction(q)
    r_vec = p.l .* dir

    # bob velocity = pivot vel + ω×r
    ωxr = cross(ω, r_vec)
    xp_dot = x_dot + ωxr[1]
    yp_dot = y_dot + ωxr[2]
    zp_dot = z_dot + ωxr[3]

    # position of bob
    x_p = x + r_vec[1]
    y_p = y + r_vec[2]
    z_p = z + r_vec[3]

    # Setpoint positions - maintaining verticality
    X = 0.0
    Y = 0.0
    Z = p.l  

    # PID error - ensure pendulum stays above pivot
    ex = x_p - (x + X)  
    ey = y_p - (y + Y)  
    ez = z_p - (z + Z)
    
    # Error velocity - control the pendulum bob's absolute motion
    dex = xp_dot 
    dey = yp_dot
    dez = zp_dot

    # ---------- Traditional PID control forces ----------
    Fx_pid = -p.Kp * ex - p.Ki * σx - p.Kd * dex
    Fy_pid = -p.Kp * ey - p.Ki * σy - p.Kd * dey
    Fz_pid = -p.Kp * ez - p.Ki * σz - p.Kd * dez
    
    # ---------- Energy-based control component ----------
    # Calculate current energy
    # Kinetic energy
    T_pivot = 0.5 * p.mpivot * (x_dot^2 + y_dot^2 + z_dot^2)
    T_bob = 0.5 * p.mp * (xp_dot^2 + yp_dot^2 + zp_dot^2)
    T_rot = 0.5 * p.mp * p.l^2 * (ω[1]^2 + ω[2]^2 + ω[3]^2)
    
    # Potential energy (reference at z=0)
    V_pivot = p.mpivot * p.g * z
    V_bob = p.mp * p.g * (z + r_vec[3])
    
    # Total current energy
    E_current = T_pivot + T_bob + T_rot + V_pivot + V_bob
    
    # Desired energy (pendulum upright, no velocity)
    V_des_pivot = p.mpivot * p.g * z
    V_des_bob = p.mp * p.g * (z + p.l)  # Bob is directly above pivot
    E_desired = V_des_pivot + V_des_bob
    
    # Energy error
    dE = E_current - E_desired
    
    # Energy control gain
    Ke = 5.0
    
    # Energy-based control forces
    F_energy_mag = -Ke * dE
    
    # Distribute energy control force based on position errors
    error_norm = sqrt(ex^2 + ey^2 + ez^2)
    if error_norm > 1e-10
        Fx_energy = F_energy_mag * (ex / error_norm)
        Fy_energy = F_energy_mag * (ey / error_norm)
        Fz_energy = F_energy_mag * (ez / error_norm)
    else
        Fx_energy = 0.0
        Fy_energy = 0.0
        Fz_energy = 0.0
    end
    
    # ---------- Combine control strategies ----------
    # Weights for blending control approaches
    w_pid = 1.0    # Weight for PID control
    w_energy = 0.0 # Weight for energy control
    
    # Combined control forces
    Fx = w_pid * Fx_pid + w_energy * Fx_energy
    Fy = w_pid * Fy_pid + w_energy * Fy_energy
    Fz = w_pid * Fz_pid + w_energy * Fz_energy
    
    # Add gravity compensation
    gravity_comp = SVector{3,Float64}(0.0, 0.0, (p.mp + p.mpivot) * p.g)
    
    # Final control forces with gravity compensation
    F = SVector{3,Float64}(Fx, Fy, Fz) + gravity_comp
    Dp = p.disturbance(t, u)  # External disturbance

    # drag
    drag_bob = p.Bp .* SVector(sign(xp_dot)*xp_dot^2,
                               sign(yp_dot)*yp_dot^2,
                               sign(zp_dot)*zp_dot^2)
    drag_piv = p.Bpivot .* SVector(sign(x_dot)*x_dot^2,
                                   sign(y_dot)*y_dot^2,
                                   sign(z_dot)*z_dot^2)

    # --- rotational dynamics ---
    I3 = I(3) # identity 3×3
    inertia_tensor = p.mp*p.l^2*(I3 - dir*dir') + ε*I3

    # gravity on bob
    gravity_force = SVector(0.0, 0.0, -p.mp * p.g)

    # net force on bob
    F_bob = gravity_force .- drag_bob .+ Dp

    # rotational dynamics 
    torque = cross(r_vec, F_bob)
    L = inertia_tensor * ω
    gyro = cross(ω, L)
    torque_damped = torque .- p.Brot .* ω
    angular_accel = inertia_tensor \ (torque_damped .- gyro)

    # translational coupling
    centripetal = cross(ω, ωxr)
    tangential = cross(angular_accel, r_vec)

    # rod‐reaction on pivot is *negative* of bob's internal acceleration
    reaction_force = -p.mp .* (centripetal .+ tangential)

    # pivot acceleration
    m_tot = p.mpivot + p.mp
    pivot_acc = (F .- drag_piv .+ reaction_force) ./ m_tot
    pivot_acc -= SVector(0.0, 0.0, p.g)

    # --- quaternion kinematics ---
    Ω = @SMatrix [
       0.0   -ω[1]  -ω[2]  -ω[3];
      ω[1]    0.0    ω[3]  -ω[2];
      ω[2]   -ω[3]   0.0    ω[1];
      ω[3]    ω[2]  -ω[1]   0.0
    ]
    q_dot = 0.5 * (Ω * q)

    # --- pack derivatives ---
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
    du[14] = ex     # integral updates
    du[15] = ey
    du[16] = ez
end

# ----------------------------- Quaternion Renormalization Callback -------------
function renormalize_q!(integrator)
    u = integrator.u
    q = SVector{4,Float64}(u[4], u[5], u[6], u[7])
    q_norm = norm(q)
    
    # Only normalize if we have a valid norm to avoid NaN issues
    if q_norm > 1e-10
        integrator.u[4:7] = q / q_norm
    end
end

# Apply callback for quaternion normalization
#cb = DiscreteCallback((u,t,integrator) -> true, renormalize_q!)

# Create a callback that normalizes approximately every 0.1 time units
time_interval = 0.00001
last_callback_time = Ref(-time_interval)  # Reference to store the last time

condition = function(u, t, integrator)
    if t - last_callback_time[] >= time_interval
        last_callback_time[] = t
        return true
    end
    return false
end
cb = DiscreteCallback(condition, renormalize_q!)

# ----------------------------- Set up & Solve ODE Problem -----------------------------
# Convert initial spherical coordinates to quaternion - start almost vertical
q0 = euler_to_quaternion(0.01, 0.01) 

# Initial condition with zero initial velocities
z0 = vcat([0.0, 0.0, 0.0], q0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
tspan = (0.0, 15.0)

# ----------------------------- Run Simulation -----------------------------
prob = ODEProblem(pendulum_quaternion!, z0, tspan, params)
sol = solve(prob, CVODE_BDF(), abstol=1e-6, reltol=1e-6, maxiters=10^6,
            callback=cb, saveat=0.05, dtmax=0.1)

println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------- Analysis Functions -----------------------------
# Extract pendulum bob position for analysis
function extract_bob_position(sol)
    positions = []
    for i in 1:length(sol.t)
        u = sol.u[i]
        x, y, z = u[1:3]
        q = normalize_quaternion(SVector{4,Float64}(u[4], u[5], u[6], u[7]))
        dir = quaternion_to_direction(q)
        r_vec = params.l .* dir
        x_p = x + r_vec[1]
        y_p = y + r_vec[2]
        z_p = z + r_vec[3]
        push!(positions, (x_p, y_p, z_p))
    end
    return positions
end

# Calculate system energy for each timestep
function calculate_energy_history(sol, params)
    energy_current = Float64[]
    energy_desired = Float64[]
    energy_error = Float64[]
    
    for i in 1:length(sol.t)
        u = sol.u[i]
        x, y, z = u[1:3]
        q = normalize_quaternion(SVector{4,Float64}(u[4], u[5], u[6], u[7]))
        x_dot, y_dot, z_dot = u[8:10]
        ω = SVector{3,Float64}(u[11], u[12], u[13])
        
        # Calculate rod direction and bob position/velocity
        dir = quaternion_to_direction(q)
        r_vec = params.l .* dir
        
        # Bob velocity
        ωxr = cross(ω, r_vec)
        xp_dot = x_dot + ωxr[1]
        yp_dot = y_dot + ωxr[2]
        zp_dot = z_dot + ωxr[3]
        
        # Calculate kinetic energy
        T_pivot = 0.5 * params.mpivot * (x_dot^2 + y_dot^2 + z_dot^2)
        T_bob = 0.5 * params.mp * (xp_dot^2 + yp_dot^2 + zp_dot^2)
        T_rot = 0.5 * params.mp * params.l^2 * (ω[1]^2 + ω[2]^2 + ω[3]^2)
        
        # Calculate potential energy
        V_pivot = params.mpivot * params.g * z
        V_bob = params.mp * params.g * (z + r_vec[3])
        
        # Current and desired energy
        E_current = T_pivot + T_bob + T_rot + V_pivot + V_bob
        E_desired = params.mpivot * params.g * z + params.mp * params.g * (z + params.l)
        
        push!(energy_current, E_current)
        push!(energy_desired, E_desired)
        push!(energy_error, E_current - E_desired)
    end
    
    return energy_current, energy_desired, energy_error
end

# Visualization function with energy analysis
function visualize_pendulum_with_energy(sol)
    # Extract positions and calculate errors
    bob_positions = extract_bob_position(sol)
    x_pos = [p[1] for p in bob_positions]
    y_pos = [p[2] for p in bob_positions]
    z_pos = [p[3] for p in bob_positions]
    
    # Calculate energy history
    energy_current, energy_desired, energy_error = calculate_energy_history(sol, params)
    
    # Create figure
    fig = Figure(resolution=(1000, 1200))  # Taller figure to accommodate additional plot
    
    # 3D trajectory plot
    ax1 = Axis3(fig[1, 1], aspect=:data, title="3D Pendulum Simulation")
    lines!(ax1, x_pos, y_pos, z_pos, color=:blue, label="Pendulum")
    
    # Add pivot trajectory
    pivot_x = [sol.u[i][1] for i in 1:length(sol.t)]
    pivot_y = [sol.u[i][2] for i in 1:length(sol.t)]
    pivot_z = [sol.u[i][3] for i in 1:length(sol.t)]
    lines!(ax1, pivot_x, pivot_y, pivot_z, color=:green, label="Pivot")
    
    # Final positions
    scatter!(ax1, [x_pos[end]], [y_pos[end]], [z_pos[end]], color=:red, markersize=20)
    scatter!(ax1, [pivot_x[end]], [pivot_y[end]], [pivot_z[end]], color=:green, markersize=10)
    
    # Position error plot
    ax2 = Axis(fig[2, 1], title="Pendulum Position Error", xlabel="Time (s)", ylabel="Error")
    t = sol.t
    error_x = [bob_positions[i][1] - sol.u[i][1] for i in 1:length(t)]
    error_y = [bob_positions[i][2] - sol.u[i][2] for i in 1:length(t)]
    error_z = [bob_positions[i][3] - (sol.u[i][3] + params.l) for i in 1:length(t)]
    
    lines!(ax2, t, error_x, color=:blue, label="X error")
    lines!(ax2, t, error_y, color=:red, label="Y error")
    lines!(ax2, t, error_z, color=:green, label="Z error")
    axislegend(ax2, position=:rt)  # Right-top position
    
    # Energy plot
    ax3 = Axis(fig[3, 1], title="System Energy", xlabel="Time (s)", ylabel="Energy (J)")
    lines!(ax3, t, energy_current, color=:blue, label="Current")
    lines!(ax3, t, energy_desired, color=:red, label="Desired")
    axislegend(ax3, position=:rt)
    
    # Energy error plot
    ax4 = Axis(fig[4, 1], title="Energy Error", xlabel="Time (s)", ylabel="Error (J)")
    lines!(ax4, t, energy_error, color=:purple)
    
    return fig
end

# Comprehensive error analysis with energy
function print_error_statistics_with_energy(sol)
    positions = extract_bob_position(sol)
    n = length(sol.t)
    
    # Calculate position errors
    error_x = Float64[]
    error_y = Float64[]
    error_z = Float64[]
    
    for i in 1:n
        x, y, z = sol.u[i][1:3]
        x_p, y_p, z_p = positions[i]
        
        push!(error_x, x_p - x)
        push!(error_y, y_p - y)
        push!(error_z, z_p - (z + params.l))
    end
    
    # Calculate energy metrics
    energy_current, energy_desired, energy_error = calculate_energy_history(sol, params)
    
    # Print position error statistics
    println("Error Statistics:")
    println("X - Max: ", maximum(abs.(error_x)), " Mean: ", sum(abs.(error_x))/n)
    println("Y - Max: ", maximum(abs.(error_y)), " Mean: ", sum(abs.(error_y))/n)
    println("Z - Max: ", maximum(abs.(error_z)), " Mean: ", sum(abs.(error_z))/n)
    
    # Print energy statistics
    println("\nEnergy Statistics:")
    println("Initial Energy: ", energy_current[1])
    println("Desired Energy: ", energy_desired[1])
    println("Final Energy: ", energy_current[end])
    println("Max Energy Error: ", maximum(abs.(energy_error)), " Mean: ", sum(abs.(energy_error))/n)
    
    # Print final position and error
    println("\nFinal Position:")
    println("Pivot: ", sol.u[end][1:3])
    println("Bob: ", positions[end])
    println("Final Position Error: (", error_x[end], ", ", error_y[end], ", ", error_z[end], ")")
    println("Final Energy Error: ", energy_error[end])
    
    # Print vertical alignment
    final_q = normalize_quaternion(SVector{4,Float64}(sol.u[end][4:7]))
    final_dir = quaternion_to_direction(final_q)
    vertical_alignment = final_dir[3]  # z-component of direction vector
    println("\nFinal Vertical Alignment: ", vertical_alignment, 
            " (1.0 = perfectly upright, -1.0 = perfectly inverted)")
end

# Run both analysis methods
print_error_statistics_with_energy(sol)

# Uncomment to visualize
fig = visualize_pendulum_with_energy(sol)
display(fig)

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
