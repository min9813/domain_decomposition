import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D


def draw(path):
    pde_data = pd.read_csv(path, sep=',')
    x1_value = pde_data['time'].astype(float)
    x2_value = [float(col[2:]) for col in pde_data.columns if 'x' in col]
    x1, x2 = np.meshgrid(x1_value, x2_value)

    x_col_name = [col for col in pde_data.columns if 'x' in col]
    z_value = pde_data[x_col_name].values

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x1, x2, z_value.T, cmap='bwr', linewidth=0)
    fig.colorbar(surf)
    ax.set_title("Surface Plot")
    plt.show()


def main(args):
    draw(args.pde_output_path)


def make_config():
    parser = argparse.ArgumentParser(
        description='PyTorch for deep face recognition')
    parser.add_argument('--pde_output_path', default='pde_output.txt')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main(make_config())
