using Plots, Measures, LinearAlgebra
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

# Function to convert from spherical to Cartesian coordinates directly
function spherical_to_cartesian(theta, phi, length)
    x = length * sin(phi) * cos(theta)
    y = length * sin(phi) * sin(theta)
    z = length * cos(phi)
    return x, y, z
end

# Function to convert Euler angles (θ, φ) to quaternion
function euler_to_quaternion(theta, phi)
    # Calculate half-angles
    theta_2 = theta/2
    phi_2 = phi/2
    
    # Calculate quaternion components
    qw = cos(phi_2) * cos(theta_2)
    qx = sin(phi_2) * cos(theta_2)
    qy = sin(phi_2) * sin(theta_2)
    qz = cos(phi_2) * sin(theta_2)
    
    return [qw, qx, qy, qz]
end

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

# Calculate pendulum position using traditional spherical coordinates
pendulum_x, pendulum_y, pendulum_z = spherical_to_cartesian(theta, phi, length)

# Calculate the projection point on the xy-plane
projection_x = pendulum_x
projection_y = pendulum_y
projection_z = 0

# Convert to quaternion (for visualization purposes only)
q = euler_to_quaternion(theta, phi)
q = q / norm(q)  # Normalize to ensure it's a unit quaternion

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

# Add quaternion visualization - showing how it relates to the same position
qw, qx, qy, qz = q
q_scale = 0.8  # Scale factor for quaternion visualization

# Add quaternion box
rect_x = -2.3
rect_y = -1.8
rect_w = 1.8
rect_h = 0.8
plot!([rect_x, rect_x+rect_w, rect_x+rect_w, rect_x, rect_x], 
      [rect_y, rect_y, rect_y+rect_h, rect_y+rect_h, rect_y],
      [0, 0, 0, 0, 0], 
      linewidth = 1, color = :black, linestyle = :dash, alpha = 0.5)

# Add quaternion text
annotate!([(rect_x+rect_w/2, rect_y+rect_h+0.1, 0, text("Quaternion Representation", 11))])
annotate!([(rect_x+0.2, rect_y+rect_h-0.2, 0, text("\$q_w = $(round(qw, digits=3))\$", 10))])
annotate!([(rect_x+0.2, rect_y+rect_h-0.4, 0, text("\$q_x = $(round(qx, digits=3))\$", 10))])
annotate!([(rect_x+0.2, rect_y+rect_h-0.6, 0, text("\$q_y = $(round(qy, digits=3))\$", 10))])
annotate!([(rect_x+rect_w-0.3, rect_y+rect_h-0.2, 0, text("\$q_z = $(round(qz, digits=3))\$", 10))])

# Draw a line from quaternion box to pendulum
connector_start_x = rect_x + rect_w/2
connector_start_y = rect_y + rect_h/2
plot!([connector_start_x, pendulum_x/2], [connector_start_y, pendulum_y/2], [0, pendulum_z/2], 
      linewidth = 1, color = :purple, linestyle = :dot)

# Display the plot
display(p)

