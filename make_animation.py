import sys
import imageio
import glob
import os
from natsort import natsorted
from structure_function import structure_function
from matplotlib import rc
rc('text', usetex=True)  # LaTeX labels

fname = sys.argv[1]
n_start = int(sys.argv[2]) if len(sys.argv) >= 2 else 0
do_delete = int(sys.argv[3]) if len(sys.argv) == 3 else 0
path = '/media/zade/Seagate Expansion Drive/Summer_Project_2019/'
filename = path + 'figs/' + fname
mhd = 'mhd' in fname or 'cgl' in fname
max_n = len(glob.glob(path+fname+'/*.athdf'))

print('Reading images from', filename)
if do_delete:
    print('Deleting images after')

if input('Do you want to continue? (Y/N) ')[0].lower() == 'n':
    print('Cancelling')
    sys.exit(0)


# Generate png images
for n in range(n_start, max_n+1):
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


# Glue together outputted image files
animate_struct()
if mhd:
    animate_mhd()

# Cleanup
if do_delete:
    for image in glob.glob(filename + '/*.png'):
        os.remove(image)
    os.rmdir(filename)
