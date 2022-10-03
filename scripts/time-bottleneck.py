import matplotlib.pyplot as plt
import numpy as np
import re
import math

INPUT_SAMPLE = 'torus-torus-torus_fallNoShift_NH_BE_interiorPoint_20220928204810t12'


def _flatten(arr):
    return np.array(np.concatenate(arr).flat)


def filepath(abs_path: str):
    """Get the path to a file relative to the current working directory.

    :param abs_path: path relative to the root directory.
    :return:
    """
    return '../' + abs_path


def output_path(fname: str):
    """Get the path to the output file relative to the working directory.

    :param fname: output file name
    :return:
    """
    return './out/' + fname


def read_file(filepath: str):
    with open(filepath) as f:
        lines = f.readlines()
        metadata = lines.pop(0), lines.pop(0), lines.pop(0)

        activities = []
        curr_activities = []
        times = []
        curr_times = []

        stat_data = []

        regex_sci_num = '\d+\.?\d*[eE]?-?\d*'
        regex_num_act = '^\d+ activities:$'
        regex_time2act = f'^\s*{regex_sci_num} s: .+$'
        regex_total = '^Total:$'
        regex_stat = '^(avg|max) # collision pairs / (iter|end of ts) = \d+$'
        for line in lines:
            if re.match(regex_num_act, line):
                # New group of activities
                if curr_activities:
                    activities.append(curr_activities)
                    times.append(curr_times)
                curr_activities = []
                curr_times = []
            elif re.match(regex_time2act, line):
                # One item in activity group
                time = float(re.findall('\d+\.?\d*[eE]?-?\d*', line)[0].replace(' s', ''))
                act = re.findall(': .+$', line)[0].replace(': ', '')
                if act == "Total":
                    continue
                curr_activities.append(act)
                curr_times.append(time)
            else:
                if curr_activities:
                    activities.append(curr_activities)
                    times.append(curr_times)
                    curr_activities = []
                    curr_times = []
                if re.match(regex_stat, line):
                    stat_data.append(line)

        return metadata, activities, times, stat_data


def plot_bars(activities, times):
    num_rows = len(activities) // 2
    num_cols = math.ceil(len(activities) // num_rows)
    fig, axs = plt.subplots(num_rows, num_cols)
    fig.suptitle(f'Bar Graphs for {INPUT_SAMPLE}')

    for i in range(len(activities)):
        y_pos = np.arange(len(activities[i]))

        subplot = axs[i // num_cols][i % num_cols]

        subplot.barh(y_pos, times[i])
        subplot.set_yticks(y_pos, labels=activities[i])

        subplot.invert_yaxis()  # labels read top-to-bottom
        subplot.set_xlabel('Time (s)')
        subplot.set_title(f'Group {i}')

    plt.show()


def plot_bars_all(act2time_tuples):
    values = [t for a, t in act2time_tuples]
    labels = [a for a, t in act2time_tuples]

    y_pos = np.arange(len(values))

    fig, ax = plt.subplots()
    # ax.barh(values, labels=labels, autopct='%.2f%%', radius=1)
    hbars = ax.barh(y_pos, values)
    ax.set_yticks(y_pos, labels=labels)

    ax.invert_yaxis()  # labels read top-to-bottom

    ax.bar_label(hbars, labels=[v if v > 0 else '' for v in values],
                 padding=8, color='b')

    ax.set_xlabel('Time (s)')
    ax.set_title(f'Bar Graph for {INPUT_SAMPLE}')
    plt.show()


def plot_pies(activities, times):
    num_rows = len(activities) // 2
    num_cols = math.ceil(len(activities) // num_rows)
    fig, axs = plt.subplots(num_rows, num_cols)
    fig.suptitle(f'Pie Graphs for {INPUT_SAMPLE}')

    for i in range(len(activities)):
        act_list = list(zip(activities[i], times[i]))
        # Data
        values = [t for a, t in act_list if t > 0]
        labels = [a for a, t in act_list if t > 0]

        subplot = axs[i // num_cols][i % num_cols]

        subplot.pie(values, labels=labels, autopct='%.2f%%', radius=1)
        subplot.set_title(f'Group {i}')

    plt.show()


def plot_pies_all(act2time_tuples):
    values = [t for a, t in act2time_tuples if t > 0]
    labels = [a for a, t in act2time_tuples if t > 0]

    fig, ax = plt.subplots()
    ax.pie(values, labels=labels, autopct='%.2f%%', radius=1)
    ax.set_title(f'Pie Graph for {INPUT_SAMPLE}')
    plt.show()


if __name__ == "__main__":
    metadata, activities, times, stat_data = read_file(
        filepath(f'output/{INPUT_SAMPLE}/info.txt'))

    # plot_pies(activities, times)
    # plot_bars(activities, times)

    flatten_times = _flatten(times)
    flatten_labels = _flatten(activities)

    act2time_tuples = list(zip(flatten_labels, flatten_times))

    plot_pies_all(act2time_tuples)
    plot_bars_all(act2time_tuples)
