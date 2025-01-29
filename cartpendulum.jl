using GLMakie
using Observables
using DifferentialEquations

# ----------------------------- System Parameters -----------------------------
const mc = 15.0      # mass of the cart
const mp = 5.0       # mass of the pendulum
const lp = 1.0       # length of the pendulum
const g = 9.81       # gravitational acceleration
const Bc = 5.0       # cart damping coefficient
const Bp = 0.5       # pendulum pivot damping coefficient

struct p
    mc::Float64  
    mp::Float64  
    lp::Float64  
    g::Float64   
    Bc::Float64  
    Bp::Float64  
end
params = p(mc, mp, lp, g, Bc, Bp)

# ----------------------------- Cart-Pendulum Model -----------------------------
function cartpend_model!(dz, z, p, t)
    # Unpack state variables
    z1, z2, z3, z4 = z

    # Equation matrices
    M = [(p.mp + p.mc)  (-p.mp*p.lp*cos(z3));
        (-p.mp*p.lp*cos(z3))  (p.mp*p.lp^2)]

    B = [(u - p.mp*p.lp*z2^2*sin(z3) - p.Bc*z2); 
        (p.mp*p.g*p.lp*sin(z3) - p.Bp*z4)]

    xdot, thetadot = (inv(M)*B)
    # Compute derivatives
    dz[1] = z2
    dz[2] = xdot
    dz[3] = z4
    dz[4] = thetadot
end

# ----------------------- Create the ODE Problem & Solve ----------------------
# Initial conditions
z0 = [0.0, 0.0, 0.5, 0.0] 
tspan = (0.0, 35.0)

# External force input
u = 0

# Solve the model ODEs
prob = ODEProblem(cartpend_model!, z0, tspan, params)
sol  = solve(prob, Rosenbrock23(), abstol=1e-8, reltol=1e-8)

# Verify the solution structure
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)

# ----------------------------- Create Figure + Plots -----------------------------
# Create cart and pendulum visualization
fig1 = Figure(size=(600, 600), fontsize=12)
ax = fig1[1, 1] = Axis(fig1, xlabel = "x", ylabel = "y", limits = (-1.5*lp, 1.5*lp, -1.5*lp, 1.5*lp), aspect=DataAspect())

cart_plot = scatter!(ax, [z0[1]], [0.0], markersize = 4*mc, color = :red) 
pole_line = lines!(ax, [z0[1], lp*sin(z0[3])], [0.0, lp*cos(z0[3])], color = :blue, linewidth = 3)
pole_plot = scatter!(ax, [lp*sin(z0[3])], [lp*cos(z0[3])], markersize = 4*mp, color = :blue)  
# display(fig1)
display(GLMakie.Screen(), fig1)

# Create phase portraits of cart and pendulum
fig2 = Figure(size=(1000, 500))
ax_pend = fig2[1, 1] = Axis(fig2, xlabel = "theta (rad)", ylabel = "thetadot (rad/s)", limits = (0,6,-7,7), title = "Pendulum Phase Portrait")
ax_cart = fig2[1, 2] = Axis(fig2, xlabel = "x (m)", ylabel = "xdot (m/s)", limits = (0,1,-1,1), title = "Cart Phase Portrait")

pendulum_line = lines!(ax_pend, Float64[], Float64[], color = :blue)
cart_line = lines!(ax_cart, Float64[], Float64[], color = :red)
#display(fig2)
display(GLMakie.Screen(), fig2)

# P1
state_obs = Observable(z0)
time_obs = Observable(0.0)
# P2
theta = Float64[]
thetadot = Float64[]
x = Float64[]
xdot = Float64[]

# ----------------------------- Update Function + visual -----------------------------
fps = 60
dt_frame = 1/fps
t_end = sol.t[end]

@async begin
    t_sim = 0.0
    
    while t_sim <= t_end
        # Sample the solution at the current simulation timestep
        z = sol(t_sim)
        
        # Update the plot
        z1, z2, z3, z4 = z
        
        # Cart locations
        cart_plot[1] = [z1]
        cart_plot[2] = [0.0]
        
        # Pendulum locations
        x_pend = [z1, z1 + lp*sin(z3)]
        y_pend = [0.0, lp*cos(z3)]
        pole_line[1] = x_pend
        pole_line[2] = y_pend
        pole_plot[1] = [z1 + lp*sin(z3)]
        pole_plot[2] = [lp*cos(z3)]

        # Store state variables
        push!(theta, z3)
        push!(thetadot, z4)
        push!(x, z1) 
        push!(xdot, z2) 
        
        # Refresh the lines on each axis
        pendulum_line[1] = theta
        pendulum_line[2] = thetadot
        cart_line[1] = x
        cart_line[2] = xdot
        
        # Update Observables
        time_obs[] = t_sim
        state_obs[] = z

        sleep(dt_frame)
        t_sim += dt_frame
    end
end

println("Continuous-time cart–pendulum simulation (ODE) is running!")


