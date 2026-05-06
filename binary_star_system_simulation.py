from vpython import *
#Web VPython 3.2

G = 6.7e-11 # universal grav constant

large_s = sphere(pos=vector(-1.1e11, 0, 0),
    mass=2e24,
    radius=7e9,
    color=color.orange,
    make_trail=True)
    
small_s = sphere(pos=vector(1.5e11, 0, 0),
    mass=1e24,
    radius=6*6.4e8,
    color=color.blue,
    make_trail=True)
    
# initializing momenta
small_s.p = small_s.mass * vector(0, 10, 0)  # velocity vector
large_s.p = -small_s.p

# individual momenta
gr_p = graph(title='Figure 1. Large/Small Stars momenta',
        xtitle=' time',
        ytitle=' p(t) ',
        align='left')
        
lspx = gcurve(graph=gr_p,
        color=color.red,
        label='large star: x-momentum')
        
lspy = gcurve(graph=gr_p,
        color=color.orange,
        label='large star: y-momentum')

lspz = gcurve(graph=gr_p,
        color=color.yellow,
        label='large star: z-momentum')
        
# small star momenta
sspx = gcurve(graph=gr_p,
        color=color.blue,
        label='small star: x-momentum')
        
sspy = gcurve(graph=gr_p,
        color=color.cyan,
        label='small star: y-momentum')

sspz = gcurve(graph=gr_p,
        color=color.green,
        label='small star: z-momentum')
        
# total momentum plot
gr_pot = graph(title='Figure 2. Total p(t) Diagram',
        xtitle='time ',
        ytitle=' p(t)',
        align='right')
        
pxtot = gcurve(graph=gr_pot,
        color=color.purple,
        label='total x-momentum in the system')
        
pytot = gcurve(graph=gr_pot,
        color=color.orange,
        label='total y-momentum in the system')
        
pztot = gcurve(graph=gr_pot,
        color=color.green,
        label='total z-momentum in the system')
        
# work and energy
gr_w = graph(title='Figure 3. Work, Kinetic Energy, Potential Energy',
        xtitle=' time ',
        ytitle=' [J] ',
        align='left')
        
w = gcurve(graph=gr_w,
    color=color.green,
    label='total internal work in system')
    
dk = gcurve(graph=gr_w,
    color=color.blue,
    label='change in kinetic energy in system')
    
du = gcurve(graph=gr_w,
    color=color.red,
    label='change in potential energy in system')
    
# total energy
gr_e = graph(title='Figure 4. Total Energy',
        xtitle=' time ',
        ytitle=' [J] ',
        align='right')
        
total = gcurve(graph=gr_e,
        color=color.green,
        label='total energy in system')
k = gcurve(graph=gr_e,
        color=color.blue,
        label='total kinetic energy in system')
u = gcurve(graph=gr_e,
        color=color.red,
        label='total potential energy in system')
        
# energy vs. separation
r_graph = graph(title=' Figure 5. Energy vs. Separation',
            xtitle=' [m] ',
            ytitle=' [J] ',
            align='left')
            
totalr = gcurve(graph=r_graph,
        color=color.green,
        label='total energy vs r')
kr = gcurve(graph=r_graph,
        color=color.blue,
        label='total kinetic energy vs r')
ur = gcurve(graph=r_graph,
        color=color.red,
        label='total potential energy vs r')
        
t = 0
dt = 30*24*60*60  # 1 month in seconds

while t < 10000 * 365 * 24 * 60 * 60:
    rate(15000)
    
    r = small_s.pos - large_s.pos
    r_mag = mag(r)

    # defining force of gravity
    Force_on_small_s = -G * large_s.mass * small_s.mass / r_mag**2 * norm(r)
    Force_on_large_s = -Force_on_small_s

    # updating momentum
    small_s.p += Force_on_small_s * dt
    large_s.p += Force_on_large_s * dt
    
    small_s.v = small_s.p / small_s.mass
    large_s.v = large_s.p / large_s.mass

    Work = (Force_on_large_s).dot(large_s.v * dt) - (Force_on_small_s).dot(small_s.v * dt)
    K = (small_s.p.mag2) / (2 * small_s.mass) + (large_s.p.mag2) / (2 * large_s.mass)
    U = -G * (small_s.mass * large_s.mass) / r_mag

    # updating position
    small_s.pos += small_s.v * dt
    large_s.pos += large_s.v * dt

    sspx.plot(pos=(t, small_s.p.x))
    sspy.plot(pos=(t, small_s.p.y))
    sspz.plot(pos=(t, small_s.p.z))

    lspx.plot(pos=(t, large_s.p.x))
    lspy.plot(pos=(t, large_s.p.y))
    lspz.plot(pos=(t, large_s.p.z))
    
    # total momentum plot
    pxtot.plot(pos=(t, small_s.p.x + large_s.p.x))
    pytot.plot(pos=(t, small_s.p.y + large_s.p.y))
    pztot.plot(pos=(t, small_s.p.z + large_s.p.z))
    
    # work plot
    w.plot(pos=(t, Work))
    dk.plot(pos=(t, K))
    du.plot(pos=(t, U))
    
    total.plot(pos=(t, K + U))
    k.plot(pos=(t, K))
    u.plot(pos=(t, U))
    
    totalr.plot(pos=(r_mag, K + U))
    kr.plot(pos=(r_mag, K))
    ur.plot(pos=(r_mag, U))

    t += dt
    
    small_s.old_v = small_s.v
    large_s.old_v = large_s.v
