import matplotlib.pyplot as plt
import numpy as np
import re
import math

INPUT_SAMPLES = {
    # 'chain5': 'torus-torus-torus_fallNoShift_NH_BE_interiorPoint_20221004145812t12',
    # 'chain10': 'torus-torus-torus_fallNoShift_NH_BE_interiorPoint_20221004145906t12',
    # 'chain35': 'torus-torus-torus_fallNoShift_NH_BE_interiorPoint_20221003103857t12',
    # 'chain100': 'torus-torus-torus_fallNoShift_NH_BE_interiorPoint_20221003232744t12',
    # 'tightFitCube1': 'ANSYS1e-3Corner_coarse_fixLowerHalf_NH_BE_interiorPoint_20221004004035t1',
    # 'tightFitCube8': 'ANSYS1e-3Corner_coarse_fixLowerHalf_NH_BE_interiorPoint_20221004003915t8',
    # 'tightFitCube12': 'ANSYS1e-3Corner_coarse_fixLowerHalf_NH_BE_interiorPoint_20221003234942t12',
    'dolphin-100': 'dolphin5K_dragright_NH_BE_interiorPoint_20230131115316t12',
    # "matOnBoard-200": 'mat40x40-mat40x40_fall_NH_BE_interiorPoint_20230127145601t12',
}


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
        # Metadata
        # 0:
        # 1: num_timesteps, num_inner_iterations,
        # 2:
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
                if act == "Total" or act == 'descent':
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


def plot_bars(input_name, activities, times):
    num_rows = len(activities) // 2
    num_cols = math.ceil(len(activities) // num_rows)
    fig, axs = plt.subplots(num_rows, num_cols)
    fig.suptitle(f'Bar Graphs for {input_name}')

    for i in range(len(activities)):
        y_pos = np.arange(len(activities[i]))

        subplot = axs[i // num_cols][i % num_cols]

        subplot.barh(y_pos, times[i])
        subplot.set_yticks(y_pos, labels=activities[i])

        subplot.invert_yaxis()  # labels read top-to-bottom
        subplot.set_xlabel('Time (s)')
        subplot.set_title(f'Group {i}')


def plot_bars_all(input_name, act2time_tuples):
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
    ax.set_title(f'Bar Graph for {input_name}')


def plot_pies(input_name, activities, times):
    num_rows = len(activities) // 2
    num_cols = math.ceil(len(activities) // num_rows)
    fig, axs = plt.subplots(num_rows, num_cols)
    fig.suptitle(f'Pie Graphs for {input_name}')

    for i in range(len(activities)):
        act_list = list(zip(activities[i], times[i]))
        # Data
        values = [t for a, t in act_list if t > 0]
        labels = [a for a, t in act_list if t > 0]

        subplot = axs[i // num_cols][i % num_cols]

        subplot.pie(values, labels=labels, autopct='%.2f%%', radius=1)
        subplot.set_title(f'Group {i}')


def plot_pies_all(input_name, act2time_tuples):
    """Plot pie graphs for one sample.

    :param input_name:
    :param act2time_tuples:
    :return:
    """
    values = [t for a, t in act2time_tuples if t > 0]
    labels = [a for a, t in act2time_tuples if t > 0]

    # Aggregate activities that is less than 5%
    total_time = sum(values)
    threshold = total_time * 0.03 or 0
    other_time = sum([t for a, t in act2time_tuples if t < threshold])
    exclude_indices = [i for i, t in enumerate(values) if t < threshold]

    def filter_excluded(ls):
        return list(map(lambda iv: iv[1], filter(lambda iv: iv[0] not in exclude_indices, enumerate(ls))))

    values = filter_excluded(values)
    labels = filter_excluded(labels)
    values.append(other_time)
    labels.append('Others')

    fig, ax = plt.subplots()
    ax.pie(values, labels=labels, autopct='%.2f%%', radius=1)
    ax.set_title(f'Pie Graph for {input_name}')


def plot_lines(activities_lst, times_lst):
    """Plot line graphs for all samples.

    :param activities_lst:
    :param times_lst:
    :return:
    """
    samples = times_lst.keys()
    sample_nums = []
    for k in samples:
        sample_nums.append(int(k.replace('chain', '')))

    # Total time
    fig_total, axs_total = plt.subplots(2)
    fig_total.suptitle('Total Running Time Plots')

    total_times = {}
    for name, times in times_lst.items():
        total_times[name] = times.sum()

    axs_total[0].set_title('Total Time vs Input Name')
    axs_total[0].plot(samples, total_times.values())
    axs_total[1].set_title('Total Time vs Number of Chains')
    axs_total[1].plot(sample_nums, total_times.values())

    # Time of each activity

    flatten_labels = list(activities_lst.values())[0]
    target_labels = ['descent',
                     'numericalFactorization',
                     'lineSearch_other',
                     'CCD',
                     'computeConstraintSets',
                     'looping EE CCS']

    # label_indices = {}
    # for label in target_labels:
    #     label_indices[label] = list(flatten_labels).index(label)

    sample2label2time = {}
    for sample in samples:
        sample2label2time[sample] = {}
        sample_acts = list(activities_lst[sample])
        for label in target_labels:
            idx = list(sample_acts).index(label)
            sample2label2time[sample][label] = times_lst[sample][idx]

    fig_acts, axs_acts = plt.subplots(len(target_labels))
    fig_acts.suptitle('Activity Running TIme Plots')

    for i in range(len(target_labels)):
        label_times = []
        label = target_labels[i]
        # for times in times_lst.values():
        #     label_times.append(times[label_indices[label]])
        for sample in samples:
            label_times.append(sample2label2time[sample][label])

        axs_acts[i].set_title(f'Time for {label} vs Number of Chains')
        axs_acts[i].plot(sample_nums, label_times)


if __name__ == "__main__":
    metadata_lst, activities_lst, times_lst, stat_data_lst = {}, {}, {}, {}
    # Plot each sample
    for sample_name, sample_file in INPUT_SAMPLES.items():
        metadata, activities, times, stat_data = read_file(
            filepath(f'output/{sample_file}/info.txt'))

        # plot_pies(activities, times)
        # plot_bars(activities, times)

        flatten_times = _flatten(times)
        flatten_labels = _flatten(activities)

        act2time_tuples = list(zip(flatten_labels, flatten_times))

        plot_pies_all(sample_name, act2time_tuples)
        plot_bars_all(sample_name, act2time_tuples)

        metadata_lst[sample_name] = metadata
        activities_lst[sample_name] = activities
        times_lst[sample_name] = times
        stat_data_lst[sample_name] = stat_data

    # Plot line graph
    # chain_acts = {}
    # chain_times = {}
    # for k in activities_lst.keys():
    #     if 'chain' in k:
    #         chain_acts[k] = _flatten(activities_lst[k])
    #         chain_times[k] = _flatten(times_lst[k])
    # plot_lines(chain_acts, chain_times)

    plt.show()
