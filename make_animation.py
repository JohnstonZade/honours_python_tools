import matplotlib.pyplot as plt
import matplotlib.animation as animation
from structure_function import structure_function
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels


def animate(n):
    print('n =', n)
    l, u, t = structure_function(fname, n)

    line1.set_data(l, u)
    line2.set_data(l, l**(2/3))
    ax.relim()
    ax.autoscale_view()
    plt.title(r"Structure Function at $t=$ " + t)
    plt.xlabel(r'log($l$)')
    plt.ylabel('log Structure Function')
    plt.legend(['Structure Function', r'$l^{2/3}$'])
    return line1, line2


fig, ax = plt.subplots()
line1, = plt.loglog([], [])
line2, = plt.loglog([], [])

fname = 'hydro_cont_turb_128'
max_n = 40

ani = animation.FuncAnimation(fig, animate, range(max_n), blit=True,
                              interval=500, repeat_delay=2000)
ani.save('animate/'+fname+'.gif')
