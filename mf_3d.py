from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np

from build.mf_3d import MagneticField

def py_parse_mf(mf_file:str):
    with open(mf_file, 'r') as file:
        content = file.read().strip().split('\n')
        content  = np.array([line.strip().split() for line in content if line[0].isdigit()]).astype(float)

    X  = content[:, 0]
    Y  = content[:, 1]
    Z  = content[:, 2]
    Bx = content[:, 3]
    By = content[:, 4]
    Bz = content[:, 5]

    fig = plt.figure()
    fig = plt.figure(figsize=(12, 45))
    ax = fig.add_subplot(projection='3d')
    ax.set_box_aspect((7.0, 1.0, 3))
    # ax.set_aspect('equalxy')
    ax.scatter(X, Z, Bz, c=Bz, cmap='viridis', s=100)
    plt.savefig('basic_pms_vis_X-Z-Bz.png')

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('file', type=str, help='mf file')

    return parser.parse_args()

if __name__ == '__main__':
    mf_file = parse_args().file

    # py_parse_mf(mf_file)

    mf = MagneticField(mf_file)

    print(mf.get_field(0.1, 0.25, 0.333))