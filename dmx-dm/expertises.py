import numpy as np
import matplotlib.pyplot as plt
import readData as rddt
import general_tools as gt
import os

def process_scanSquids(dirPath, fb0, pixel_start, pixel_end):
    xmin, xmax = -2**10, 2**10
    ymin, ymax = -2**13, 2**13
    ymin3, ymax3 = -0.12, 0.12
    xlabel = 'Feedback (ADU)'
    ylabel1 = 'Error (ADU)'
    ylabel2 = 'System gain (ADU/ADU)'
    ylabel3 = 'kNorm'

    plotDirname = "PLOTS"
    pathPlot = os.path.join(dirPath, plotDirname)
    gt.createdir(pathPlot)
    plotFileName = os.path.join(pathPlot, 'scanSquids.png')

    files = [f for f in os.listdir(dirPath) \
             if os.path.isfile(os.path.join(dirPath, f)) \
             and f[-5:] == ".fits"]

    if len(files) == 0:
        raise ValueError("No FITS files.")
    if len(files) > 1:
        raise ValueError("Too much FITS files.")

    print("Processing file ", files[0])
    xName, ctrl, feedback, error = rddt.readScan(os.path.join(dirPath, files[0]))
    error = np.mean(error[pixel_start:pixel_end, :], axis=0); # keeping the data of a single pixel

    if xName != 'Feedback':
        raise ValueError("Wrong fits file type.")

    print("Found {0:} steps in the scan".format(len(error)))

    # computing the local gain (derivative of the function)
    derivative = gt.derivative(feedback, error)

    # Looking for operating point data
    diff = np.abs(feedback - fb0)
    ifb0 = np.where(diff == diff.min())
    error0 = error[ifb0].mean()
    fb0 = feedback[ifb0].mean()
    gain0 = derivative[ifb0].mean()

    # Doing the plot
    fig = plt.figure(figsize=(10, 13))
    fig.suptitle('Scan Feedback ' + files[0].split('_')[1], fontsize=16)
    loc = 'upper right'

    ax0 = fig.add_subplot(4, 1, 1)
    mksz, mksz2 = 2, 7
    ax0.plot(feedback, error, '.b', markersize=mksz)
    ax0.set_ylim(ymin, ymax)
    ax0.grid()
    ax0.set_xlabel(xlabel)
    ax0.set_ylabel(ylabel1)

    ax1 = fig.add_subplot(4, 1, 2)
    mksz, mksz2 = 2, 7
    ax1.plot(feedback, error, '.b', markersize=mksz)
    text = "fb0={0:}, Lockpoint={1:}".format(int(fb0), int(error0))
    ax1.plot(fb0, error0, 'or', markersize=mksz2, label=text)
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax1.grid()
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1)
    ax1.legend(loc=loc)

    ax2 = fig.add_subplot(4, 1, 3)
    ax2.plot(feedback[:-1], derivative, '.g', markersize=mksz)
    text = "fb0={0:}, gain={1:6.2f}".format(int(fb0), gain0)
    ax2.plot(fb0, gain0, 'or', markersize=mksz2, label=text)
    ax2.set_xlim(xmin, xmax)
    ax2.grid()
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel2)
    ax2.legend(loc=loc)

    ax3 = fig.add_subplot(4, 1, 4)
    ax3.plot(feedback[:-1], 1/derivative, '.k', markersize=mksz)
    text = "fb0={0:}, knorm={1:6.3f}".format(int(fb0), 1/gain0)
    ax3.plot(fb0, 1/gain0, 'or', markersize=mksz2, label=text)
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(ymin3, ymax3)
    ax3.grid()
    ax3.set_xlabel(xlabel)
    ax3.set_ylabel(ylabel3)
    ax3.legend(loc=loc)

    for item in ([ax0.xaxis.label, ax0.yaxis.label, ax1.xaxis.label, ax1.yaxis.label,
                  ax2.xaxis.label, ax2.yaxis.label, ax3.xaxis.label, ax3.yaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(10)
    for item in (ax0.get_xticklabels() + ax0.get_yticklabels() + ax1.get_xticklabels() + ax1.get_yticklabels()
                 + ax2.get_xticklabels() + ax2.get_yticklabels() + ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(10)

    fig.tight_layout()
    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')


def plot_scan(dirPath, pixel):

    plotDirname = "PLOTS"
    pathPlot = os.path.join(dirPath, plotDirname)
    gt.createdir(pathPlot)

    files = [f for f in os.listdir(dirPath) \
             if os.path.isfile(os.path.join(dirPath, f)) \
             and f[-5:] == ".fits"]

    if len(files) == 0:
        raise ValueError("No FITS files.")
    if len(files) > 1:
        raise ValueError("Too much FITS files.")

    print("Processing file ", files[0])
    xName, ctrl, xscan, error = rddt.readScan(os.path.join(dirPath, files[0]))
    error = error[pixel,:] # keeping the data of a single pixel
    plotFileName = os.path.join(pathPlot, 'scan_'+xName+'.png')

    # Doing the plot
    fig = plt.figure(figsize=(10, 6))
    fig.suptitle('Scan ' + xName + ' ' + files[0].split('_')[1], fontsize=16)
    loc = 'upper right'

    ax1 = fig.add_subplot(1, 1, 1)
    mksz = 2
    ax1.plot(xscan, error, '.b', markersize=mksz)
    ax1.grid()
    ax1.set_xlabel(xName + ' (ADU)')
    ax1.set_ylabel('Error (ADU)')

    for item in ([ax1.xaxis.label, ax1.yaxis.label]):
        item.set_weight('bold')
        item.set_fontsize(10)
    for item in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(10)

    fig.tight_layout()
    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')


dirPath = ['/Users/laurent/Data/TestPlan20-TDM/scanSquids20241206_120125',
    '/Users/laurent/Data/TestPlan20-TDM/scan_20250130-105328_Feedback'
]
pixel = 0
fb0 = -528
for path in dirPath:
    process_scanSquids(path, fb0, 1,33)

dirPath = [
    '/Users/laurent/Data/TestPlan20-TDM/scan_20250129-162135_Offset',
    '/Users/laurent/Data/TestPlan20-TDM/scan_20250129-162403_Feedback',
    '/Users/laurent/Data/TestPlan20-TDM/scan_20250130-104725_Feedback',
    '/Users/laurent/Data/TestPlan20-TDM/scan_20250130-104901_Feedback',
    '/Users/laurent/Data/TestPlan20-TDM/scan_20250130-104949_Feedback',
    '/Users/laurent/Data/TestPlan20-TDM/scan_20250130-105328_Feedback'
]
pixel = 0
for path in dirPath:
    plot_scan(path, pixel)
