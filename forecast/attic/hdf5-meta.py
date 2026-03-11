import h5py

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filename", "-f", help="HDF5 file", required=True)
args = parser.parse_args()

f = h5py.File(args.filename, 'r')

for k in f.attrs.keys():
    print(f"{k} => {f.attrs[k]}")
