"""
Call with the command line arguments
    [1] sim -- SIMBA or Illustris
    [2] suite -- LH or 1P
    [3] redshift index (e.g. 004)
"""

from sys import argv
import os.path
import subprocess
import pathlib
from time import time

SIM = argv[1]
SUITE = argv[2]
ZIDX = argv[3]
RECOMPUTE = bool(int(argv[4]))

assert SIM in ['SIMBA', 'Illustris']
assert SUITE in ['LH', '1P', 'CV']
assert int(ZIDX) >= 0 and int(ZIDX) <= 33

NSIMS = 1000 if SUITE=='LH' \
        else 66 if SUITE=='1P' \
        else 27 if SUITE=='CV' \
        else None

INPUT_ROOT = '/projects/QUIJOTE/CAMELS/Sims/%s/%s'%('IllustrisTNG' if SIM=='Illustris' else 'SIMBA',
                                                    SUITE)
OUTPUT_ROOT = '/projects/QUIJOTE/Leander/CAMELS_meanT/%s/%s'%('IllustrisTNG' if SIM=='Illustris' else 'SIMBA',
                                                                 SUITE)

FNAMES = ['grp_M200c.bin', 'grp_meanT.bin', ]

for idx in range(NSIMS) :
    
    input_dir = '%s_%d'%(INPUT_ROOT, idx)
    fgrp = os.path.join(input_dir, 'fof_subhalo_tab_%s.hdf5'%ZIDX)
    fprt = os.path.join(input_dir, 'snap_%s.hdf5'%ZIDX)

    assert os.path.isfile(fgrp)
    assert os.path.isfile(fprt)

    output_dir = '%s_%d/%s'%(OUTPUT_ROOT, idx, ZIDX)

    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

    # check if this one has already been done
    if not RECOMPUTE and any(os.path.isfile(os.path.join(output_dir, f)) for f in FNAMES) :
        assert all(os.path.isfile(os.path.join(output_dir, f)) for f in FNAMES)
        print('%s already done'%output_dir)
        continue

    # otherwise we compute
    tstart = time()

    try :
        subprocess.run(['./examples/CAMELS_profiles_%s'%SIM.upper(), fgrp, fprt, output_dir],
                       check=True)
    except subprocess.CalledProcessError :
        print('%s excited with non-zero exit code. Continuing.'%output_dir)

    tend = time()
    print('Did %s in %f seconds.'%(output_dir, tend-tstart))

    assert all(os.path.isfile(os.path.join(output_dir, f)) for f in FNAMES)
