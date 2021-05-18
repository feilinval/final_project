  
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt	
from matplotlib.animation import FuncAnimation
	
frames = 365
seconds_in_year = 365 * 24 * 60 * 60
years = 1
t = np.linspace(0, years*seconds_in_year, frames)

def move_func(s, t):
    (x, v_x, y, v_y,
     xm, vm_x, ym, vm_y,
     xq, vq_x, yq, vq_y,
     xf, vf_x, yf, vf_y,
     xh, vh_x, yh, vh_y) = s
    
    dxdt = v_x
    dv_xdt = - G * m * x / (x**2 + y**2)**1.5
    dydt = v_y
    dv_ydt = - G * m * y / (x**2 + y**2)**1.5
    
    dxmdt = vm_x
    dvm_xdt = - G * m * xm / (xm**2 + ym**2)**1.5
    dymdt = vm_y
    dvm_ydt = - G * m * ym / (xm**2 + ym**2)**1.5
    
    dxqdt = vq_x
    dvq_xdt = - G * m * xq / (xq**2 + yq**2)**1.5
    dyqdt = vq_y
    dvq_ydt = - G * m * yq / (xq**2 + yq**2)**1.5

    dxfdt = vf_x
    dvf_xdt = - G * m * xf / (xf**2 + yf**2)**1.5
    dyfdt = vf_y
    dvf_ydt = - G * m * yf / (xf**2 + yf**2)**1.5
    
    dxhdt = vh_x
    dvh_xdt = - G * m * xh / (xh**2 + yh**2)**1.5
    dyhdt = vh_y
    dvh_ydt = - G * m * yh / (xh**2 + yh**2)**1.5

    return (dxdt, dv_xdt, dydt, dv_ydt,
            dxmdt, dvm_xdt, dymdt, dvm_ydt,
            dxqdt, dvq_xdt, dyqdt, dvq_ydt,
            dxfdt, dvf_xdt, dyfdt, dvf_ydt,
            dxhdt, dvh_xdt, dyhdt, dvh_ydt
            )

G = 6.67 * 10**(-11)
m = 1.98 * 10**(30)
ae = 149 * 10**9

x0 = ae
v_x0 = 0
y0 = 0
v_y0 = 30000

xm0 = 1.523662 * ae
vm_x0 = 0
ym0 = 0
vm_y0 = 24000

xq0 = 0.38709927 * ae
vq_x0 = 0
yq0 = 0
vq_y0 = 47360

xf0 = 0.723332 * ae
vf_x0 = 0
yf0 = 0
vf_y0 = 35020

xh0 = 5.204267 * ae
vh_x0 = 0
yh0 = 0
vh_y0 = 13070

s0 = (x0, v_x0, y0, v_y0, 
      xm0, vm_x0, ym0, vm_y0,
      xq0, vq_x0, yq0, vq_y0,
      xf0, vf_x0, yf0, vf_y0,
      xh0, vh_x0, yh0, vh_y0)

def solve_func(i, key):
    sol = odeint(move_func, s0, t)
    if key == 'point':
        x = sol[i, 0]
        y = sol[i, 2]
        xm = sol[i,4]
        ym = sol[i,6]
        xq = sol[i, 8]
        yq = sol[i, 10]
        xf = sol[i, 12]
        yf = sol[i, 14]
        xh = sol[i, 16]
        yh = sol[i, 18]
    else:
        x = sol[:i, 0]
        y = sol[:i, 2]
        xm =sol[:i, 4]
        ym = sol[:i, 6]
        xq = sol[:i, 8]
        yq = sol[:i, 10]
        xf = sol[:i, 12]
        yf = sol[:i, 14]
        xh = sol[:i, 16]
        yh = sol[:i, 18]
    return ((x, y), (xm, ym), (xq, yq), (xf, yf),(xh,yh))

fig, ax = plt.subplots()

ball, = plt.plot([], [], 'o', color='b')
ball_line, = plt.plot([], [], '-', color='b')

ballm, = plt.plot([], [], 'o', color='b')
ball_linem, = plt.plot([], [], '-', color='b')

ballq, = plt.plot([], [], 'o', color='b')
ball_lineq, = plt.plot([], [], '-', color='b')

ballf, = plt.plot([], [], 'o', color='b')
ball_linef, = plt.plot([], [], '-', color='b')

ballh, = plt.plot([], [], 'o', color='b')
ball_lineh, = plt.plot([], [], '-', color='b')

plt.plot([0], [0], 'o', color='y', ms=15)

def animate(i):
    ball.set_data(solve_func(i, 'point')[0])
    ball_line.set_data(solve_func(i, 'line')[0])

    ballm.set_data(solve_func(i, 'point')[1])
    ball_linem.set_data(solve_func(i, 'line')[1])
    
    ballq.set_data(solve_func(i, 'point')[2])
    ball_lineq.set_data(solve_func(i, 'line')[2])

    ballf.set_data(solve_func(i, 'point')[3])
    ball_linef.set_data(solve_func(i, 'line')[3])

    ballh.set_data(solve_func(i, 'point')[4])
    ball_lineh.set_data(solve_func(i, 'line')[4])

ani = FuncAnimation(fig,
                    animate,
                    frames=frames,
                    interval=30)  

plt.axis('equal')
edge = 6 * ae
ax.set_xlim(-edge, edge)
ax.set_ylim(-edge, edge)

plt.show()
