from analyse import analyse

fnames = ['cgl_31Jan_320160_nuc1', 'cgl_31Jan_320160_nuc4', 'cgl_31Jan_320160_nuc10', 'cgl_28Jan_24048_nuc1e4']

for fname in fnames:
    print('Running', fname)
    analyse(fname)
