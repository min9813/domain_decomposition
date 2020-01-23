import argparse
import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.animation as animation
import yaml
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

DRAW_Y_VAL = [0.1, 0.3, 0.5, 0.7, 0.9]


def load_result(args):
    pde_data_list = []
    params = {}
    use_exp_name = []
    for exp_v in args.demo_exp:
        config_path = "./config/{}.yaml".format(exp_v)
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        # if config["dim"] != 1:
            # raise ValueError
        if len(params) == 0:
            params = config
        elif args.same_shape:
            assert params["delta"]["t"] == config["delta"]["t"], exp_v
            assert params["delta"]["x"] == config["delta"]["x"], exp_v

            assert np.all(params["time"] == config["time"]), exp_v
            assert np.all(params["space"]["x"] == config["space"]["x"]), exp_v

        path = "./result/{}/{}".format(exp_v, exp_v)
        print("load from {}".format(path))
        pde_data = pd.read_csv(path, sep=',')
        pde_data_list.append(pde_data)
        if args.base_version in exp_v:
            true_idx = args.demo_exp.index(exp_v)
        use_exp_name.append(exp_v.split("_")[0][3:])

    use_exp_name = '_'.join(use_exp_name)
    return pde_data_list, use_exp_name, true_idx


def draw_1d(path):
    pde_data = pd.read_csv(path, sep=',')
    x1_value = pde_data['time'].astype(float)
    x2_value = [float(col[2:]) for col in pde_data.columns if 'x' in col]
    x1, x2 = np.meshgrid(x1_value, x2_value)

    x_col_name = [col for col in pde_data.columns if 'x' in col]
    z_value = pde_data[x_col_name].values

    print(z_value.shape)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x1, x2, z_value.T, cmap='bwr', linewidth=0)
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    ax.set_zlabel('u')

    fig.colorbar(surf)
    # name_dict = {
    # "exp5": "exact",
    # "exp6"
    # }
    # ax.set_title("2d result")
    exp_dir = os.path.dirname(path)
    img_save_path = os.path.join(exp_dir, "result_all.png")
    plt.savefig(img_save_path)


def calc_1d_error(args):
    pde_data_list, use_exp_name, true_idx = load_result(args)
    true_df = pde_data_list[true_idx]
    x_cols = [col for col in true_df.columns if 'x' in col]
    true_values = true_df[x_cols].values
    mask = true_df["time"].isin(pde_data_list[-1]["time"])
    print(true_df["time"])
    use_time = true_df["time"][mask]
    true_values = true_values[mask]

    for idx, df in enumerate(pde_data_list):
        if idx == true_idx:
            continue
        mean_arr = np.mean(np.abs(true_values - df[x_cols].values))
        print(args.labels[idx], ":", mean_arr)

    if args.demo_exp[0] == "exp5":
        true_values = pde_data_list[1][x_cols].values
        mask = pde_data_list[1]["time"].isin(pde_data_list[-1]["time"])
        for idx, df in enumerate(pde_data_list):
            if idx in (true_idx, 1):
                continue
            mean_arr = np.mean(np.abs(true_values - df[x_cols].values))
            print(args.labels[idx], ":", mean_arr)


def draw_1d_video(args):
    pde_data_list, use_exp_name, true_idx = load_result(args)
    pde_data = pde_data_list[0]
    t_value = pde_data['time'].astype(float)
    x_cols = [col for col in pde_data.columns if 'x' in col]
    x_value = [float(col[2:]) for col in x_cols]

    result_dir = './result/1d/'
    if os.path.exists(result_dir) is False:
        os.makedirs(result_dir)
    save_demo_path = os.path.join(
        result_dir, 'all_ani_{}_{}.mp4'.format(use_exp_name, args.diff_method))
    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    # turn to h,w
    # video_shape = tuple(video_shape.tolist())
    last_index = np.inf

    fig = plt.figure()
    shape_array = (fig.bbox_inches._points[1]*100).astype(int).tolist()
    # print('shape:', shape_array)
    # video_shape = tuple(shape_array)
    video_shape = (1600, 1200)
    plt.close()
    print('save to', save_demo_path)
    video = cv2.VideoWriter(save_demo_path, fourcc, 10., video_shape)
    colors = ["green", "skyblue", "orange", "black"]
    NUM = 4
    for idx, now_time in tqdm(enumerate(t_value), total=len(t_value)):
        fig = plt.figure(figsize=(8.0, 6.0))
        w, h = (fig.bbox_inches._points[1]*100).astype(int).tolist()
        ax1 = fig.add_subplot(111)
        diff_fig = [plt.figure(figsize=(8.0, 6.0)) for _ in range(NUM-1)]
        diff_ax = [f.add_subplot(111) for f in diff_fig]
        exact_value = pde_data_list[true_idx].loc[pde_data["time"] ==
                                                  now_time, x_cols].values.astype(float).reshape(-1)
        for _idx, pde_data in enumerate(pde_data_list):
            field_value = pde_data.loc[pde_data["time"] ==
                                       now_time, x_cols].values.astype(float).reshape(-1)
            if args.diff_method == "raw":
                mask = np.ones_like(exact_value[1:-1]).astype(bool)
                diff = field_value - exact_value
                title = "sim - true "
                y_lim = (-0.15, 0.1)
            elif args.diff_method == "raw_relative":
                mask = exact_value != 0
                diff = np.zeros_like(exact_value)
                diff[mask] = (field_value - exact_value)[mask] / \
                    exact_value[mask]
                title = "sim - true rel "
                y_lim = (-0.5, 0.5)
            else:
                raise ValueError

            ax1.plot(x_value, field_value,
                     label=args.labels[_idx], color=colors[_idx])
            try:
                # print(diff)
                diff_ax[_idx].bar(x_value, diff, width=0.05,
                                  color=colors[_idx], align="center")
                diff_ax[_idx].scatter(
                    x_value, diff, marker='o', linewidths="0.5", color=colors[_idx])
                # diff_ax[_idx].text(0.35, 0.4, "error={:.4f}".format(
                #     np.mean(np.abs(diff[mask]))), size=20, color="black")
                diff_ax[_idx].set_title(
                    "diff {} {}".format(title, args.labels[_idx]))
                diff_ax[_idx].set_xlabel("x")
                # diff_ax[_idx].set_ylabel("diff")
                diff_ax[_idx].set_xlim(-0.1, 1.1)
                diff_ax[_idx].set_ylim(*y_lim)
                # diff_ax[_idx].legend()
            except IndexError:
                continue

        ax1.set_title("result t={:.4f}".format(now_time))
        ax1.set_xlabel("x")
        ax1.set_ylabel("u")
        ax1.set_xlim(-0.1, 1.1)
        ax1.set_ylim(-0.1, 1.1)
        ax1.legend()
        per_num = np.sqrt(NUM).astype(int)
        canvas = np.zeros((video_shape[1], video_shape[0], 3)).astype(np.uint8)
        for i in range(per_num):
            for j in range(per_num):
                canvas_idx = i * per_num + j
                if canvas_idx == 0:
                    buf, size = fig.canvas.print_to_buffer()

                else:
                    buf, size = diff_fig[canvas_idx-1].canvas.print_to_buffer()

                img_arr = np.frombuffer(buf, dtype=np.uint8).reshape(
                    size[1], size[0], -1)
                im = cv2.cvtColor(img_arr, cv2.COLOR_RGBA2BGR)
                canvas[h*i:h*(i+1), w*j:w*(j+1)] = im
                if canvas_idx > 0:
                    diff_fig[canvas_idx-1].clear()
        video.write(canvas)
        plt.close()
        fig.clear()

        if idx > last_index - 1:
            break

    video.release()


def calc_2d_error(args):
    PREFIX = "(x y)="
    pde_data_list, use_exp_name, true_idx = load_result(args)
    true_df = pde_data_list[true_idx]
    x_cols = [col for col in true_df.columns if PREFIX in col]
    true_values = true_df[x_cols].values
    use_time = true_df["time"]
    # print(true_df["time"])
    # print(pde_data_list[-1]["time"])
    # print(pde_data_list[1]["time"])

    for idx, df in enumerate(pde_data_list):
        mask1 = use_time.isin(df["time"])
        mask2 = df["time"].isin(use_time)
        if idx == true_idx:
            continue
        mean_arr = np.mean(
            np.abs(true_values[mask1] - df[x_cols].values[mask2]))
        print(args.labels[idx], ":", mean_arr)

    if args.demo_exp[0] == "exp1":
        true_values = pde_data_list[0][x_cols].values
        mask = pde_data_list[1]["time"].isin(pde_data_list[-1]["time"])
        for idx, df in enumerate(pde_data_list):
            if idx in (true_idx, 0):
                continue
            mean_arr = np.mean(np.abs(true_values - df[x_cols].values))
            print(args.labels[idx], ":", mean_arr)


def draw_2d(path, args):
    PREFIX = "(x y)="
    result_dir = os.path.dirname(path)

    pde_data = pd.read_csv(path, sep=',')

    xy_col = [col for col in pde_data.columns if PREFIX in col]
    xy_bar = [tuple(map(float, col.replace(PREFIX, "").split()))
              for col in xy_col]
    x_bar, y_bar = map(lambda x: np.array(x), zip(*xy_bar))
    x_col_num = len(np.unique(x_bar))
    y_col_num = len(np.unique(y_bar))
    print(np.unique(x_bar))
    print(f"surface shape: ({x_col_num}, {y_col_num})")
    x_bar = x_bar.reshape(y_col_num, x_col_num)
    y_bar = y_bar.reshape(y_col_num, x_col_num)

    last_index = np.inf
    bname = os.path.basename(result_dir)
    if os.path.exists(result_dir) is False:
        os.makedirs(result_dir)
    skip_time = max(1, len(pde_data.index) // args.video_length)
    save_demo_path = os.path.join(result_dir, f'animate_skip{skip_time}.mp4')
    # w,h of figure
    # fig = plt.figure()
    # video_shape = (fig.bbox_inches._points[1] * 100).astype(int)
    # print()

    y_diff = (np.max(y_bar) - np.min(y_bar)) / len(y_bar)
    x_diff = (np.max(x_bar) - np.min(x_bar)) / len(x_bar)

    y_lim = np.min(y_bar)-y_diff, np.max(y_bar)
    x_lim = np.min(x_bar)-x_diff, np.max(x_bar)
    u_max = np.max(pde_data.loc[0, xy_col])
    u_min = 0 - u_max * 0.1

    x_diff = x_lim[1] - x_lim[0]
    ratio = x_diff / 6
    x_size = 10
    y_size = 5

    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    # turn to h,w
    # video_shape = tuple(video_shape.tolist())
    print("fig size:", (x_size, y_size))
    fig = plt.figure(figsize=(x_size, y_size))
    shape_array = (fig.bbox_inches._points[1]*100).astype(int).tolist()
    video_shape = tuple(shape_array)
    plt.close()
    print('save to', save_demo_path)
    video = cv2.VideoWriter(save_demo_path, fourcc, 10., video_shape)

    if skip_time > 1:
        max_log = np.log(np.max(pde_data.index)+1)
        min_log = np.log(np.min(pde_data.index)+1)
        max_log = np.arange(min_log, max_log+1, (max_log -
                                                 min_log) / args.video_length)[:args.video_length]
        index = np.unique(np.exp(max_log).astype(int))
        mask = index <= np.max(pde_data.index)
        index = index[mask]
        index = np.concatenate((np.zeros(1), index))
    else:
        index = pde_data.index
    print(index)

    for idx in tqdm(index):
        now_time = pde_data.loc[idx, "time"]
        source_value = pde_data.loc[idx, xy_col].values.reshape(x_bar.shape)
        # print(source_value)
        # print(x_bar)
        # print(y_bar)
        # print(x_bar[0])
        # print(y_bar[0])
        # raise ValueError
        fig = plt.figure(figsize=(x_size, y_size))
        ax = fig.add_subplot(111, projection="3d")
        surf = ax.plot_surface(
            x_bar, y_bar, source_value, cmap='bwr', linewidth=0)
        plt.colorbar(surf)
        # ax.set_title("3d result {} t={:.6f}".format(bname, now_time))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("u")
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)
        ax.set_zlim(u_min, u_max)
        ax.view_init(elev=6, azim=-110)
        buf, size = fig.canvas.print_to_buffer()
        img_arr = np.frombuffer(buf, dtype=np.uint8).reshape(
            size[1], size[0], -1)
        im = cv2.cvtColor(img_arr, cv2.COLOR_RGBA2BGR)
        video.write(im)

        plt.close()

        if idx > last_index - 1:
            break

    video.release()


def concat_video(paths):
    result_dir = './result/2d/'
    if os.path.exists(result_dir) is False:
        os.makedirs(result_dir)

    cap_files = [cv2.VideoCapture(path) for path in paths]
    use_exp_name = "_".join(
        [os.path.basename(os.path.dirname(path)) for path in paths])
    for c in cap_files:
        if c.isOpened() is False:
            raise ValueError

    save_demo_path = os.path.join(
        result_dir, 'all_ani_{}_concat.mp4'.format(use_exp_name))

    w = cap_files[0].get(cv2.CAP_PROP_FRAME_WIDTH)
    h = cap_files[0].get(cv2.CAP_PROP_FRAME_WIDTH)
    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    last_index = np.inf

    video_shape = (1600, 1200)
    video = cv2.VideoWriter(save_demo_path, fourcc, 10., video_shape)

    frame_num = int(cap_files[0].get(cv2.CAP_PROP_FRAME_COUNT))

    for frame_now in tqdm(range(frame_num)):
        images = []
        for cap_file in cap_files:
            flag, image = cap_file.read()
            if flag is False:
                raise ValueError
            image = cv2.resize(image, (800, 600))
            images.append(image)

        new_img = np.concatenate(images[:2], axis=1)
        new_img2 = np.concatenate(images[2:], axis=1)
        new_img = np.concatenate((new_img, new_img2), axis=0)
        video.write(new_img)
        if frame_now > last_index:
            break

    video.release()


def draw_2d_demo(args):
    PREFIX = "(x y)="

    pde_data_list, use_exp_name, true_idx = load_result(args)
    pde_data = pde_data_list[0]

    t_value = pde_data['time'].astype(float)
    xy_col = [col for col in pde_data.columns if PREFIX in col]
    xy_bar = [tuple(map(float, col.replace(PREFIX, "").split()))
              for col in xy_col]

    use_cols = [[xy_col[idx] for y_val in DRAW_Y_VAL if col[1] == y_val]
                for idx, col in enumerate(xy_bar)]

    result_dir = './result/2d/'
    if os.path.exists(result_dir) is False:
        os.makedirs(result_dir)
    save_demo_path = os.path.join(
        result_dir, 'all_ani_{}_{}.mp4'.format(use_exp_name, args.diff_method))
    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    # turn to h,w
    # video_shape = tuple(video_shape.tolist())
    last_index = np.inf

    fig = plt.figure()
    shape_array = (fig.bbox_inches._points[1]*100).astype(int).tolist()
    # print('shape:', shape_array)
    video_shape = tuple(shape_array)
    plt.close()
    print('save to', save_demo_path)
    video = cv2.VideoWriter(save_demo_path, fourcc, 10., video_shape)
    colors = ["green", "skyblue", "orange", "black"]
    for idx, now_time in tqdm(enumerate(t_value), total=len(t_value)):
        fig = plt.figure(figsize=(8.0, 6.0))
        w, h = (fig.bbox_inches._points[1]*100).astype(int).tolist()
        ax1 = fig.add_subplot(111)

        for _idx, pde_data in enumerate(pde_data_list):
            field_value = pde_data.loc[pde_data["time"] ==
                                       now_time, x_cols].values.astype(float).reshape(-1)

            ax.plot(x_value, field_value,
                    label=args.labels[_idx], color=colors[_idx])
        ax.set_title("result t={:.4f}".format(now_time))
        ax.set_xlabel("x")
        ax.set_ylabel("u")
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.legend()
        buf, size = fig.canvas.print_to_buffer()
        img_arr = np.frombuffer(buf, dtype=np.uint8).reshape(
            size[1], size[0], -1)
        im = cv2.cvtColor(img_arr, cv2.COLOR_RGBA2BGR)
        video.write(im)
        plt.close()
        if idx > last_index - 1:
            break

    video.release()


def main(args):
    if args.mode == "1d_demo":
        if args.base_version == "exp5":
            args.demo_exp = ["exp6", "exp7", "exp8", "exp5"]
            args.labels = ["exp exp5",
                           "imp exp7", "exp_imp exp11", "exact exp9"]
        elif args.exp_version == "exp21":
            args.demo_exp = ["exp21", "exp22", "exp23"]
            args.labels = ["exact exp21", "imp exp22", "domain exp23"]
        else:
            raise ValueError
        draw_1d_video(args)
    elif args.mode == "1d_calc":
        if args.base_version == "exp5":
            args.demo_exp = ["exp5", "exp6", "exp7", "exp8", "exp9"]
            args.labels = ["exact exp5", "exp exp6",
                           "imp exp7", "exp_imp exp8", "decompose imp exp9"]
        elif args.base_version == "exp21":
            args.demo_exp = ["exp21", "exp22", "exp23"]
            args.labels = ["exact exp21", "imp exp22", "domain exp23"]
        else:
            raise ValueError
        calc_1d_error(args)
    elif args.mode == "1d_error_ana":
        if args.base_version == "exp5":
            args.demo_exp = ["exp5", "exp30", "exp31", "exp32",
                             "exp33", "exp34", "exp35", "exp36", "exp37", "exp38"]
            args.labels = ["exact exp5", "exp_imp 1", "exp_imp 2", "exp_imp 3", "exp_imp 4", "exp_imp 5", "exp_imp 6",
                           "exp_imp 7", "exp_imp 8", "exp_imp 9"]
        elif args.base_version == "exp49":
            args.demo_exp = ["exp49", "exp40", "exp41", "exp42",
                             "exp43", "exp44", "exp45", "exp46", "exp47", "exp48"]
            args.labels = ["exact exp49", "exp_imp 1.0", "exp_imp 2.0", "exp_imp 3", "exp_imp 4", "exp_imp 5", "exp_imp 6",
                           "exp_imp 7", "exp_imp 8", "exp_imp 9"]

        else:
            raise ValueError
        calc_1d_error(args)
    elif args.mode == "2d_error_ana":
        if args.base_version == "exp10_exact":
            args.demo_exp = ["exp10_exact", "exp51", "exp52", "exp53",
                             "exp54", "exp55", "exp56", "exp57", "exp58"]
            args.labels = ["exact exp49", "exp_imp 2d 1.0", "exp_imp 2.0", "exp_imp 3", "exp_imp 4", "exp_imp 5", "exp_imp 6",
                           "exp_imp 7", "exp_imp 8"]
        else:
            raise ValueError
        calc_2d_error(args)

    elif args.mode == "2d_demo":
        args.demo_exp = ["exp1", "exp2", "exp3", "exp10_exact", "exp4"]
        args.labels = ["exp exp1", "imp exp2", "exp_imp exp3",
                       "exact exp10", "decompose imp exp4"]
    elif args.mode == "2d_calc":
        if args.base_version == "exp10_exact":
            args.demo_exp = ["exp1", "exp2", "exp3", "exp10_exact", "exp4"]
            args.labels = ["exp exp1", "imp exp2",
                           "exp_imp exp3", "exact exp10", "decompose imp exp4"]

        elif args.base_version == "exp12":
            args.demo_exp = ["exp12", "exp13", "exp14"]
            args.labels = ["exp exp12", "imp exp13", "exp_imp exp14"]
        elif args.base_version == "exp15":
            args.demo_exp = ["exp15", "exp16", "exp17"]
            args.labels = ["exp exp15", "imp exp16", "exp_imp exp17"]
        elif args.base_version == "exp17":
            args.demo_exp = ["exp17", "exp18", "exp19"]
            args.labels = ["exp exp17", "imp exp18", "exp_imp exp19"]
        else:
            raise NotImplementedError
        calc_2d_error(args)

    elif args.mode == "one":
        yaml_path = "./config/{}.yaml".format(args.exp_version)
        with open(str(yaml_path), "r") as f:
            config = yaml.safe_load(f)
        args.space_dim = int(config["dim"])
        result_dir = args.exp_version
        args.pde_output_path = "./result/{}/{}".format(result_dir, result_dir)
        print("load from {}".format(args.pde_output_path))
        if args.space_dim == 1:
            draw_1d(args.pde_output_path)
        elif args.space_dim == 2:
            draw_2d(args.pde_output_path, args)
        else:
            raise NotImplementedError
    elif args.mode == "concat":
        demo_exp = ["exp1", "exp2", "exp3", "exp10_exact"]
        path = ["./result/{}/animate.mp4".format(d) for d in demo_exp]
        concat_video(path)
    else:
        raise NotImplementedError


def str2bool(ttt):
    return "t" in ttt


def make_config():
    parser = argparse.ArgumentParser(
        description='PyTorch for deep face recognition')
    parser.add_argument('--exp_version', default='exp1')
    parser.add_argument('--mode', default='demo',
                        choices=["1d_demo", "2d_demo", "one", "concat", "1d_calc", "2d_calc", "1d_error_ana", "2d_error_ana"])
    parser.add_argument('--diff_method', default='raw',
                        choices=["raw", "raw_relative"])
    parser.add_argument('--same_shape', default='t', type=str2bool)
    parser.add_argument('--base_version', default='exp5')
    parser.add_argument('--video_length', default=100, type=int)

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main(make_config())
