using Plots, Measures, LinearAlgebra, LaTeXStrings
gr() # Use GR backend for 3D plotting

# --------------------------- Spherical Formulation ---------------------------
# Set up the 3D plot with a zoomed-in view
p = plot(
    size = (1000, 800),
    xlim = (-2.5, 2.5), 
    ylim = (-2.5, 2.5),
    zlim = (-1, 2.5),
    aspect_ratio = :equal,
    legend = false,
    camera = (30, 20),
    grid = false,    
    ticks = false,      
    framestyle = :none,
    margin = -100mm  
)

# Pendulum parameters
phi = π/6    # polar angle from z-axis
theta = π/4  # azimuthal angle in xy-plane
length = 2.0 # pendulum length

# Calculate pendulum position in 3D
pendulum_x = length * sin(phi) * cos(theta)
pendulum_y = length * sin(phi) * sin(theta)
pendulum_z = length * cos(phi)

# Calculate the projection point on the xy-plane
projection_x = pendulum_x
projection_y = pendulum_y
projection_z = 0

# Draw coordinate axes with dashed lines and labels
plot!([-3, 3], [0, 0], [0, 0], linestyle = :dash, linewidth = 2, color = :red)
plot!([0, 0], [-3, 3], [0, 0], linestyle = :dash, linewidth = 2, color = :green)
plot!([0, 0], [0, 0], [-3, 2.5], linestyle = :dash, linewidth = 2, color = :blue)

# Add axis labels at the end of each axis line
annotate!([2.7], [0], [0], ["\$x\$"])
annotate!([0], [2.7], [0], ["\$y\$"])
annotate!([0], [0], [2.6], ["\$z\$"])

# Draw the pendulum rod in 3D
plot!([0, pendulum_x], [0, pendulum_y], [0, pendulum_z], linewidth = 4, color = :black, label = "Rod")

# Draw the projection lines (vertical and horizontal)
plot!([pendulum_x, projection_x], [pendulum_y, projection_y], [pendulum_z, projection_z], 
      linestyle = :dash, linewidth = 1, color = :black, alpha = 0.5)
plot!([0, projection_x], [0, projection_y], [0, projection_z], 
      linestyle = :dash, linewidth = 1, color = :black, alpha = 0.5)

# Add the masses as 3D points
scatter!([0], [0], [0], markersize = 15, color = :blue, label = "Pivot Mass")
scatter!([pendulum_x], [pendulum_y], [pendulum_z], markersize = 12, color = :red, label = "Pendulum Mass")

# Add text annotations in 3D
annotate!([(0.65, -0.4, -0.05, text("\$m_{pivot}\$", 13))])
annotate!([(pendulum_x + 0.2, pendulum_y + 0.3, pendulum_z, text("\$m_{p}\$", 13))])

# Draw and label the phi angle arc (in xz-plane)
phi_arc_radius = 0.65
phi_arc_points = 30
phi_arc = range(0, phi, length = phi_arc_points)
phi_arc_x = phi_arc_radius * sin.(phi_arc)
phi_arc_y = zeros(phi_arc_points)
phi_arc_z = phi_arc_radius * cos.(phi_arc)
plot!(phi_arc_x, phi_arc_y, phi_arc_z, linewidth = 2, color = :purple)
# Add phi label at the middle of the arc
phi_label_x = phi_arc_radius * sin(phi/2) * 1.2
phi_label_z = phi_arc_radius * cos(phi/2) * 1.2
annotate!([(phi_label_x, 0, phi_label_z, text("\$φ\$", 12))])

# Draw and label the theta angle arc (in xy-plane)
theta_arc_radius = 0.65
theta_arc_points = 30
theta_arc = range(0, theta, length = theta_arc_points)
theta_arc_x = theta_arc_radius * cos.(theta_arc)
theta_arc_y = theta_arc_radius * sin.(theta_arc)
theta_arc_z = zeros(theta_arc_points)
plot!(theta_arc_x, theta_arc_y, theta_arc_z, linewidth = 2, color = :orange)
# Add theta label at the middle of the arc
theta_label_x = theta_arc_radius * cos(theta/2) * 1.2
theta_label_y = theta_arc_radius * sin(theta/2) * 1.2
annotate!([(theta_label_x, theta_label_y, 0, text("\$θ\$", 12))])

# Display the plot
display(p)

# -------------------------- Quarternion Formulation --------------------------
u  = normalize([-0.6, 0.25, 0.7])  # rotation axis
alpha = π  # rotation angle
q0 = cos(alpha/2)  # scalar part
qv = sin(alpha/2)*u  # vector part

# rotation helper
rotate(p) = (q0^2 - dot(qv,qv)).*p .+
            2*(dot(qv,p)).*qv .+
            2*q0*cross(qv,p)

# Spherical Pendulum Parameters
φ = π/4 * 180/π # polar angle
θ = π/4 * 180/π # azimuthal
L = 1.0  # length

# original pendulum tip & projection
tip = [L*sin(φ)*cos(θ), L*sin(φ)*sin(θ), L*cos(φ)]
proj = [tip[1], tip[2], 0.0]
tip_r  = rotate(tip)
proj_r = rotate(proj)
pivot  = [0.0,0.0,0.0]

# Build 3D Figure
p = plot(
  size        = (1000,800),
  xlim        = (-1.5,1.5), ylim = (-1.5,1.5),
  zlim        = (-0.5,1.5),
  aspect_ratio= :equal,
  legend      = false,
  camera      = (30,20),
  grid        = false,
  ticks       = false,
  framestyle  = :none,
  margin      = -100mm
)

# — Dashed coordinate axes —
plot!([-3,3],[0,0],[0,0], linestyle=:dash, lw=2, color=:red)
plot!([0,0],[-3,3],[0,0], linestyle=:dash, lw=2, color=:green)
plot!([0,0],[0,0],[-3,2.0], linestyle=:dash, lw=2, color=:blue)
annotate!([(1.7,0,0,text("\$x\$",12)),
           (0,1.7,0,text("\$y\$",12)),
           (0,0,1.6,text("\$z\$",12))])

# Original pendulum (solid)
plot!([0, tip[1]], [0, tip[2]], [0, tip[3]], lw=4, color=:black)
scatter!([pivot[1]], [pivot[2]], [pivot[3]], markersize=15, color=:blue)
scatter!([tip[1]], [tip[2]], [tip[3]], markersize=12, color=:red)
annotate!([(0.55, -0.3, -0.0125, text("\$m_{pivot}\$",13)),(tip[1]+0.15, tip[2]+0.15, tip[3], text("\$m_p\$",13))])

# Rotated pendulum (dashed + lighter)
# rod
plot!([0, tip_r[1]], [0, tip_r[2]], [0, tip_r[3]], lw=3, color=:black, linestyle=:dash, alpha=0.4)
scatter!([pivot[1]], [pivot[2]], [pivot[3]], markersize=15, color=:blue, alpha=0.3)
scatter!([tip_r[1]], [tip_r[2]], [tip_r[3]], markersize=12, color=:red, alpha=0.3)

# Rotation axis + angle arc + quaternion text
plot!([-u[1], u[1]], [-u[2], u[2]], [-u[3], u[3]], lw=3, color=:blue)
annotate!([(u[1]*1.1, u[2]*1.1, u[3]*1.1, text("\$Rotation  Axis\$",10,:blue))])

# arc
dir0 = tip .- (dot(tip,u))*u
w1 = normalize(dir0)            
w2 = cross(u, w1)               
arc_radius = 0.85
arc_offset = 0.2
npts = 101              
ts = range(0, stop=alpha, length=npts)
arc_pts = [ arc_offset*u .+ arc_radius*(cos(t)*w1 .+ sin(t)*w2) for t in ts]
arc = hcat(arc_pts...)
plot!(arc[1,:], arc[2,:], arc[3,:], lw = 3, color = :purple, clip = false)
mid_idx = (npts+1) ÷ 2
mid_pt  = arc_pts[mid_idx]
lbl_off = 0.1*(cos(alpha/2)*w1 .+ sin(alpha/2)*w2)
lbl_pos = mid_pt .+ lbl_off
annotate!([(lbl_pos[1], lbl_pos[2], lbl_pos[3]-0.05,text(L"\alpha", 14, :purple))])
# rotation axis label
rotation = L"\mathbf{u} = (u_1, u_2, u_3)"
annotate!([(-1.65, -1.65, 1.45, text(rotation, 12, halign = :left))])
# quaternion label
qstr = L"q = (q_0,q_1,q_2,q_3) = (\cos(\alpha/2),\;\sin(\alpha/2)\,\mathbf{u})"
annotate!([(-1.65, -1.65, 1.3, text(qstr, 12, halign = :left))])

# Display the plot
display(p)
