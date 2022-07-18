import sys
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# TODO: Include introgressed segment areas
def compare_3(naive_gamma, bw_gamma, true_intro_states):
    ng = naive_gamma[:, 1]
    bw = bw_gamma[:, 1]
    tis = true_intro_states
    x_axis_start = 0
    x_axis_end = 40

    # Zoom in on introgressed states
    # zoom_start = 8500
    # zoom_stop = 8620
    # ng = naive_gamma[zoom_start:zoom_stop, 1]
    # bw = bw_gamma[zoom_start:zoom_stop, 1]
    # tis = true_intro_states[zoom_start:zoom_stop]
    # x_axis_start = zoom_start
    # x_axis_end = zoom_stop

    # ORIGINAL THREE PLOTS
    # plt.rcParams.update({'font.size': 12})
    # fig, axes = plt.subplots(1, 3, figsize=(12, 9), sharex=False, sharey=False, dpi=100)
    # axes[0].plot(np.linspace(x_axis_start, x_axis_end, ng.shape[0]), ng, marker='.', color='tab:blue')
    # axes[0].spines['top'].set_visible(False)
    # axes[0].spines['right'].set_visible(False)
    # axes[0].set_title('Naive HMM')
    # axes[1].plot(np.linspace(x_axis_start, x_axis_end, bw.shape[0]), bw, marker='.', color='tab:green')
    # axes[1].spines['top'].set_visible(False)
    # axes[1].spines['right'].set_visible(False)
    # axes[1].set_title('Baum-Welch HMM')
    # axes[2].plot(np.linspace(x_axis_start, x_axis_end, tis.shape[0]), tis, marker='.', color='tab:orange')
    # axes[2].spines['top'].set_visible(False)
    # axes[2].spines['right'].set_visible(False)
    # axes[2].set_title('True introgressed states')
    # # Delete "Measured in thousands" when zooming
    # fig.supxlabel('500bp window measured in thousands', size=14, weight='bold')
    # fig.supylabel('Percent of introgression in window according to model', size=14, weight='bold')
    # fig.suptitle('Identifying Introgressed Segments', weight='bold')
    # fig.tight_layout()
    # plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)

    # OVERLAYING INTROGRESSED STATES - 2 PLOTS
    # plt.rcParams.update({'font.size': 12})
    # fig, axes = plt.subplots(1, 2, figsize=(12, 9), sharex=False, sharey=False, dpi=100)
    # # axes[0].plot(np.linspace(x_axis_start, x_axis_end, tis.shape[0]), tis, marker='.', color='tab:orange')
    # axes[0].plot(np.linspace(x_axis_start, x_axis_end, ng.shape[0]), ng, marker='.', color='tab:blue')
    # axes[0].spines['top'].set_visible(False)
    # axes[0].spines['right'].set_visible(False)
    # axes[0].set_title('Naive HMM')
    #
    # # axes[1].plot(np.linspace(x_axis_start, x_axis_end, tis.shape[0]), tis, marker='.', color='tab:orange')
    # axes[1].plot(np.linspace(x_axis_start, x_axis_end, bw.shape[0]), bw, marker='.', color='tab:green')
    # axes[1].spines['top'].set_visible(False)
    # axes[1].spines['right'].set_visible(False)
    # axes[1].set_title('Baum-Welch HMM')
    #
    # # Delete "Measured in thousands" when zooming
    # fig.supxlabel('500bp window measured in thousands', size=14, weight='bold')
    # fig.supylabel('Percent of introgression in window according to model', size=14, weight='bold')
    # fig.suptitle('Identifying Introgressed Segments', weight='bold')
    # fig.tight_layout()
    # plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)

    # OVERLAYING INTROGRESSED STATES: 1 PLOT BW HMM ONLY
    x_axis_start = 1
    x_axis_end = 5
    tis = np.array((0, 0, 1, 1, 0))
    plt.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(12, 9), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(np.linspace(x_axis_start, x_axis_end, tis.shape[0]), tis, marker='.', color='tab:orange')
    ax.plot(np.linspace(x_axis_start, x_axis_end, ng.shape[0]), ng, marker='.', color='tab:blue')
    ax.plot(np.linspace(x_axis_start, x_axis_end, bw.shape[0]), bw, marker='.', color='tab:green')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Observed Sequence: NNCCN')

    fig.supxlabel('Character in Sequence', size=14, weight='bold')
    fig.supylabel('HMM Confidence Locus is Introgressed', size=14, weight='bold')
    fig.suptitle('Naive HMM vs Baum-Welch HMM', weight='bold')
    fig.tight_layout()
    plt.legend(["Hidden States", "Naive HMM", "Baum-Welch HMM"], loc="upper left")
    plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)

    # # OVERLAYING INTROGRESSED STATES: 1 PLOT NAIVE HMM ONLY
    # x_axis_start = 1
    # x_axis_end = 5
    # tis = np.array((0, 0, 1, 1, 0))
    #
    # plt.rcParams.update({'font.size': 12})
    # fig = plt.figure(figsize=(12, 9), dpi=100)
    # ax = fig.add_subplot(111)
    # ax.plot(np.linspace(x_axis_start, x_axis_end, tis.shape[0]), tis, marker='.', color='tab:orange')
    # ax.plot(np.linspace(x_axis_start, x_axis_end, ng.shape[0]), ng, marker='.', color='tab:green')
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # ax.set_title("Naive HMM on Observed Sequence NNCCN")
    #
    # fig.supxlabel('Character in Sequence', size=14, weight='bold')
    # fig.supylabel('Confidence that hidden state is introgressed according to Naive HMM', size=14, weight='bold')
    # fig.suptitle('Identifying Introgressed Segments', weight='bold')
    # fig.tight_layout()
    # plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)

    # ZOOMING
    # plt.rcParams.update({'font.size': 12})
    # fig, axes = plt.subplots(1, 2, figsize=(12, 3), sharex=False, sharey=False, dpi=100)
    # axes[0].plot(np.linspace(x_axis_start, x_axis_end, bw.shape[0]), bw, marker='.', color='tab:green')
    # axes[0].spines['top'].set_visible(False)
    # axes[0].spines['right'].set_visible(False)
    # axes[0].set_title('Baum-Welch HMM')
    # axes[1].plot(np.linspace(x_axis_start, x_axis_end, ng.shape[0]), ng, marker='.', color='tab:blue')
    # axes[1].spines['top'].set_visible(False)
    # axes[1].spines['right'].set_visible(False)
    # axes[1].set_title('Naive HMM')
    # # Delete "Measured in thousands" when zooming
    # fig.supxlabel('500bp window', size=14, weight='bold')
    # fig.supylabel('Percent of introgression in window according to model', size=14, weight='bold')
    # fig.suptitle('Identifying Introgressed Segments', weight='bold')
    # fig.tight_layout()
    # plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)

# TODO: Include introgressed segment areas
# Input: performances: an numpy array of OL numpy arrays where the OLth array represents the HMM performance values
#        on this sequence at optimization limit OL
def display_performance(performances):
    ol = performances.shape[0]
    print(ol)

    # create data
    false_pos_r = performances[:, 0].T
    miss_rate = performances[:, 1].T
    sensitivity = performances[:, 2].T
    specificity = performances[:, 3].T

    # SUBPLOTTED PERFORMANCE
    # plt.rcParams.update({'font.size': 12})
    # fig, axes = plt.subplots(1, 4, figsize=(12, 9), sharex=False, sharey=False, dpi=100)
    # axes[0].plot(np.linspace(0, ol, false_pos_r.shape[0]), false_pos_r, marker='.', color='tab:red')
    # axes[0].spines['top'].set_visible(False)
    # axes[0].spines['right'].set_visible(False)
    # axes[0].set_title('False Positive Rate')
    # axes[1].plot(np.linspace(0, ol, miss_rate.shape[0]), miss_rate, marker='.', color='tab:green')
    # axes[1].spines['top'].set_visible(False)
    # axes[1].spines['right'].set_visible(False)
    # axes[1].set_title('Miss Rate')
    # axes[2].plot(np.linspace(0, ol, sensitivity.shape[0]), sensitivity, marker='.', color='tab:blue')
    # axes[2].spines['top'].set_visible(False)
    # axes[2].spines['right'].set_visible(False)
    # axes[2].set_title('Sensitivity')
    # axes[3].plot(np.linspace(0, ol, specificity.shape[0]), specificity, marker='.', color='tab:orange')
    # axes[3].spines['top'].set_visible(False)
    # axes[3].spines['right'].set_visible(False)
    # axes[3].set_title('Specificity')
    # # Delete "Measured in thousands" when zooming
    # fig.supxlabel('Number of optimization steps in Baum-Welch', size=14, weight='bold')
    # fig.supylabel('Percentage', size=14, weight='bold')
    # fig.suptitle('Measuring performance improvement using Baum-Welch with normalization and a 90% threshold', weight='bold')
    # fig.tight_layout()
    # plt.savefig('/Users/briankirz/Downloads/performance_test.png', facecolor='white', dpi=300)

    # TODO: SUBPLOTS
    # # Y-axis values
    y1 = false_pos_r
    y2 = miss_rate
    y3 = sensitivity
    y4 = specificity
    #
    # # Function to plot
    # plt.plot(y2)
    # plt.plot(y4)
    # # plt.plot(y3)
    # # plt.plot(y4)
    #
    # # Function add a legend
    # plt.legend(["Miss Rate", "Specificity"], loc="lower left") # "Sensitivity", "Specificity"], loc="lower right")
    # plt.xlabel("Number of Baum-Welch iterations")
    # plt.ylabel("Percentage")
    # plt.show()
    #
    # # SECOND SUBPLOT
    # plt.plot(y1)
    # plt.plot(y3)
    # plt.legend(["False Positive Rate", "Sensitivity"], loc="upper left")
    # plt.xlabel("Number of Baum-Welch iterations")
    # plt.ylabel("Percentage")
    # plt.show()

    fig, ax = plt.subplots(2, 2)

    ax[0, 0].plot(range(y1), 'r')  # row=0, col=0
    ax[1, 0].plot(range(y2), 'b')  # row=1, col=0
    ax[0, 1].plot(range(y3), 'g')  # row=0, col=1
    ax[1, 1].plot(range(y4), 'k')  # row=1, col=1
    plt.legend(["False Positive Rate", "Miss Rate", "Sensitivity", "Specificity"], loc="lower left")
    plt.xlabel("Number of Baum-Welch iterations")
    plt.ylabel("Percentage")
    plt.show()