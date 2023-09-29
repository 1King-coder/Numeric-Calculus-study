import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from zeros_of_functions.Functions import Func
import sympy as sym
from math import sin, cos
"""

figure, axe = plt.subplots()
axe.set(xlim=[-10, 20], ylim=[-0.5, 50], xlabel='X', ylabel='Y')

x = np.arange(-10, 10, 0.1)

f = (x - 2)*(x - 3)

func = axe.plot(x, f, label="f(x)")

f_line = [axe.plot(x, (2*i - 5)*(x - i) + (i-2)*(i-3), c='r', label="f'(x)") for i in x]
f_scatter = [axe.plot(i, (i - 2)*(i - 3), c='r',marker='.') for i in x]

"""
"""
def update_frame (frame) -> tuple:
    x_axe = x[:frame]
    y = g[:frame]
    
    tg_line.set_xdata(x_axe)
    tg_line.set_ydata(y)

    

    return [tg_line]

anime = animation.FuncAnimation(figure, update_frame, 1000, interval=10, blit=True)
"""
"""
derivative_line_animation = animation.ArtistAnimation(fig=figure, artists=f_line, interval=0)
derivative_dot_animation = animation.ArtistAnimation(fig=figure, artists=f_scatter, interval=0)
"""

x = sym.symbols('x')
func = Func((x - 2)*(x - 3), x)

res, iter_x, iter_y = func.newton_method(0, 10**(-6)).values()

iter_x, iter_y = [float(value) for value in iter_x], [float(value) for value in iter_y]

x_np = np.arange(-10, 10, 0.1)
print (type(x_np))

fig, axe = plt.subplots()
axe.set(xlim= [-1, 4], ylim= [-0.25, 1])

axe.plot(x_np, func.f(x_np))
axe.plot(x_np, [0 for  _ in x_np], c='black', label='x')
axe.plot([0 for _ in x_np], x_np, c='black', label='y')

derivative = [axe.plot(x_np, (x_np - i)*func.f(i, func.dfunc) + func.f(i), c='r') for i in iter_x]
where_derivative_touches_x = [axe.plot(i, 0, c='orange', marker='.') for i in iter_x[1:]]

axe.plot(res, 0, c='r', marker='.')
anim = animation.ArtistAnimation(fig, derivative, interval=0)

anim_2 = animation.ArtistAnimation(fig, where_derivative_touches_x, interval=0)
plt.show()