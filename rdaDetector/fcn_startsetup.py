import os

try:
    import fooof
    print ('fooof already installed...')
except ImportError:
    print ('Installing fooof...')
    os.system("pip install fooof")

try:
    import hdf5storage
    print ('hdf5storage already installed...')
except ImportError:
    print ('Installing hdf5storage...')
    os.system("pip install hdf5storage")

try:
    import mne
    print ('mne already installed...')
except ImportError:
    print ('Installing mne...')
    os.system("pip install mne")
