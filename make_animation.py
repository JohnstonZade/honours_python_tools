import sys
import imageio
import glob
import os
from natsort import natsorted
from structure_function import structure_function
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels

fname = str(sys.argv[1])
mhd = 'mhd' in fname or 'cgl' in fname
max_n = 200
filename = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
filename += 'figs/' + fname

for n in range(max_n+1):
    print('n =', n)
    structure_function(fname, n, do_mhd=mhd)


def animate_struct():
    struct = natsorted(glob.glob(filename + '/struct*.png'))
    struct_ims = []

    for image in struct:
        struct_ims.append(imageio.imread(image))

    imageio.mimsave('animate/' + fname + '_struct.gif', struct_ims,
                    duration=0.2)


def animate_mhd():
    parl = natsorted(glob.glob(filename + '/*0.png'))
    # Figure out way to get end index depending on number of angles
    perp = natsorted(glob.glob(filename + '/*2.png'))
    parl_ims, perp_ims = [], []

    for image in parl:
        parl_ims.append(imageio.imread(image))
    for image in perp:
        perp_ims.append(imageio.imread(image))

    imageio.mimsave('animate/' + fname + '_parl.gif', parl_ims,
                    duration=0.2)
    imageio.mimsave('animate/' + fname + '_perp.gif', perp_ims,
                    duration=0.2)


animate_struct()
if mhd:
    animate_mhd()

for image in glob.glob(filename + '/*.png'):
    os.remove(image)
