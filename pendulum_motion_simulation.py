from vpython import *
#Web VPython 3.2

scene = canvas(background = color.white)

fixed_point = sphere(pos = vec(0, 8, 0),
        radius = 0.2,
        color = color.orange)
        
# constants
l_0 = 8 # length of string
theta_0 = 10 * pi/180 # angle in degrees
g = 9.81 # grav constant 

# objects
ball = sphere(pos = fixed_point.pos + vec(l_0 * sin(theta_0), -l_0 * cos(theta_0), 0),
    mass = 2,
    radius = 0.5,
    color = color.red)
    
string = helix(pos = fixed_point.pos, # not gonna act as a string
    axis = ball.pos - fixed_point.pos,
    radius = 0.2,
    color = color.yellow)
    
# initializing angular variables as attributes of the ball pendulum object
ball.theta = theta_0
ball.omega = 0

# initializing momentum
ball.p = ball.mass * l_0 * ball.omega # using v = omega * length of string

# graphinhg
th_graph = graph(title = 'theta-omega', xtitle = 't', ytitle = 'y')
theta_g = gcurve(graph = th_graph, color = color.green)
omega_g = gcurve(graph = th_graph, color = color.blue)

# initializing time constant and step
t = 0
dt = 0.001

# defining period calculation variables
ball.old_omega = ball.omega
ball.T = 0 # period timer

while t < 20:
    rate(3000)
    
    F_tan = -ball.mass * g * sin(ball.theta)

    ball.p += F_tan * dt * l_0
    
    ball.omega = ball.p / (ball.mass * l_0)
    
    ball.theta += ball.omega * dt
    
    ball.pos = fixed_point.pos + vec(l_0 * sin(ball.theta), -l_0 * cos(ball.theta), 0)

    string.axis = ball.pos - fixed_point.pos

    # plotting (refer to previous lavs)
    theta_g.plot(t, ball.theta)
    omega_g.plot(t, ball.omega)
    
    # update time
    t = t + dt
    
    # update period timer
    ball.T = ball.T + dt
    
    # measure period by printing T when pendulum reaches original position
    # hint: use angualr velocity to figure this out
    
    if ball.old_omega < 0 and ball.omega > 0:
        print('Measured Period: ', ball.T, ' s')
        print('Theoretical Period: ', 2*pi*sqrt(l_0/g), ' s')
        ball.T = 0  # Reset the period counter
        
    ball.old_omega = ball.omega