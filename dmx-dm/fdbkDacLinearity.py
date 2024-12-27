# imports
import numpy as np
import matplotlib.pyplot as plt
import os
import general_tools as gt
from matplotlib.ticker import MultipleLocator
import glob
import pandas

path = "/Users/laurent/Data/TestPlan21-perfo/20241120_fdbkLinearity/"
dacNbBits = 14
nbCols = 4
xtit = "Feedback value (ADU)"
ytit1 = "Feedback value (V)"
ytit2 = "Non Linearity (% FSR)"
ylim = [-1, 1]
xlim = [-2**(dacNbBits-1), 2**(dacNbBits-1)]


def setGrid(ax, majorOn, MinorOn, xmajor, xminor, ymajor, yminor):
    # Définir le pas de la grille majeure
    ax.xaxis.set_major_locator(MultipleLocator(xmajor))
    ax.yaxis.set_major_locator(MultipleLocator(ymajor))

    # Définir le pas de la grille mineure
    ax.xaxis.set_minor_locator(MultipleLocator(xminor))
    ax.yaxis.set_minor_locator(MultipleLocator(yminor))

    # Activation de la grille majeure
    ax.grid(majorOn, which='major', linestyle='-', linewidth='0.8', color='black')  # Grille majeure

    # Activer la grille mineure
    if (MinorOn):
        ax.minorticks_on()  # Activer les ticks mineurs
        ax.grid(True, which='minor', linestyle=':', linewidth='0.5', color='gray')

# Checking the indexes are associated by pairs
def checkPairs(indx):
    diff = indx[1:] - indx[:-1]
    jumps = np.where(abs(diff) > 1)[0]
    toBeDeleted = []
    if indx[0]%2==1:
        toBeDeleted.append(0)
    for i in range(len(jumps)):
        if indx[jumps[i]+1] % 2 == 1:
            toBeDeleted.append(jumps[i]+1)
        if indx[jumps[i]] % 2 == 0:
            toBeDeleted.append(jumps[i])
    return np.delete(indx, toBeDeleted)

#a=np.array([0,1,2,3,4,5,9,10,11,12,13,16,17])
#print(checkPairs(a))

def plotLinearity(adu, meas, indx, axx):
    meas = meas[indx]
    measUp = meas[::2]
    measDown = meas[1::2]
    meas = np.concatenate((measUp, measDown))

    adu = adu[indx]
    aduUp = adu[::2]
    aduDown = adu[1::2]
    adu = np.concatenate((aduUp, aduDown))

    # Linear fit of the data
    coeffs = np.polyfit(adu, meas, 1)
    fit = coeffs[1] + coeffs[0] * adu
    fsrV = fit.max() - fit.min()
    fitUp = coeffs[1] + coeffs[0] * aduUp
    fitDown = coeffs[1] + coeffs[0] * aduDown
    deviationPercentFSRUp = 100 * (measUp - fitUp) / fsrV
    deviationPercentFSRDown = 100 * (measDown - fitDown) / fsrV

    axx.scatter(aduUp, deviationPercentFSRUp, marker='.', color='green', label='Increasing values')
    axx.scatter(aduDown, deviationPercentFSRDown, marker='.', color='blue', label='Decreasing values')
    axx.set_xlim(xlim)
    axx.set_xlabel(xtit)
    axx.set_ylabel(ytit2)
    ylim = axx.get_ylim()
    delta = ylim[1] - ylim[0]
    axx.set_ylim(ylim[0]-delta/2, ylim[1]+delta/2)
    ytit = ylim[1]+delta * 0.6

    return (min(deviationPercentFSRUp.min(), deviationPercentFSRDown.min()),
            max(deviationPercentFSRUp.max(), deviationPercentFSRDown.max()), ytit)

files = glob.glob(path + "*.csv")

if len(files) == 0:
    print("Error no csv files")
elif len(files) > 0:
    fullFileName = os.path.join(path, files[0])
    if len(files) > 1:
        print("Error wrong number of csv files")
print("Processing " + fullFileName)

plotDirname = "PLOTS"
pathPlot = os.path.join(path, plotDirname)
gt.createdir(pathPlot)

df = pandas.read_csv(fullFileName)
df.columns = ["col0", "col1", "col2", "col3"]

columns = [ np.array(df.loc[:, "col0"].astype(float)),
        np.array(df.loc[:, "col1"].astype(float)),
        np.array(df.loc[:, "col2"].astype(float)),
        np.array(df.loc[:, "col3"].astype(float))]

adu_base = np.zeros(2**dacNbBits)
for i in range(2**(dacNbBits - 1)):
    adu_base[i * 2] = -2 ** (dacNbBits - 1) + i * 2
    adu_base[i * 2 + 1] = 2 ** (dacNbBits - 1) - 1 - i * 2

for colId in range(nbCols):

    col = columns[colId]

    derivative = np.concatenate((col[:-2] - col[2:], [0,0]))
    indexes = np.where((col < 2) & (np.abs(derivative) < 0.1))[0]
    indexes = checkPairs(indexes)
    k2 = 0.03
    indexes2 = indexes[(indexes >= k2 * (2 ** (dacNbBits))) & (indexes <= (1-k2) * (2 ** (dacNbBits)))]
    indexes2 = checkPairs(indexes2)
    k3 = 0.2
    indexes3 = indexes[(indexes >= k3 * (2 ** (dacNbBits))) & (indexes <= (1-k3) * (2 ** (dacNbBits)))]
    indexes3 = checkPairs(indexes3)

    plotFileName = os.path.join(pathPlot, 'fdbkLinearity_col{0:}.png'.format(colId))

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle('Feedback linearity measurement for column {0:}'.format(colId), fontsize=16)

    ax1 = fig.add_subplot(3, 2, 2)
    min1, max1, ytit1 = plotLinearity(adu_base, col, indexes, ax1)
    ax1.text(-2**(dacNbBits-1)+256, ytit1, "Linearity on the full range: {0:3.2f} / {1:3.2f} % of FSR".format(min1, max1))

    ax2 = fig.add_subplot(3, 2, 4)
    min2, max2, ytit2 = plotLinearity(adu_base, col, indexes2, ax2)
    ax2.text(-2**(dacNbBits-1)+256, ytit2, "Linearity on 95% of the range: {0:3.2f} / {1:3.2f} % of FSR".format(min2, max2))

    ax3 = fig.add_subplot(3, 2, 6)
    min3, max3, ytit3 = plotLinearity(adu_base, col, indexes3, ax3)
    ax3.text(-2**(dacNbBits-1)+256, ytit3, "Linearity on 60% of the range: {0:3.2f} / {1:3.2f} % of FSR".format(min3, max3))

    meas = col[indexes]
    measUp = meas[::2]
    measDown = meas[1::2]
    meas = np.concatenate((measUp, measDown))

    adu = adu_base[indexes]
    aduUp = adu[::2]
    aduDown = adu[1::2]
    adu = np.concatenate((aduUp, aduDown))

    # Linear fit of the data
    coeffs = np.polyfit(adu, meas, 1)
    fit = coeffs[1] + coeffs[0] * adu

    ax0 = fig.add_subplot(1, 2, 1)
    ax0.scatter(aduUp, measUp, marker='.', color='green', label='Increasing values')
    ax0.scatter(aduDown, measDown, marker='.', color='blue', label='Decreasing values')
    ax0.plot(adu, fit, color='red', label='Linear fit')
    #ax0.set_xlim(xlim)
    #ax0.set_ylim(ylim)
    ax0.set_xlabel(xtit)
    ax0.set_ylabel(ytit1)
    ax0.legend(loc='best')
    #setGrid(ax0, True, True, 2**12, 2**8, 1, 0.2)

    fig.tight_layout()

    plt.savefig(plotFileName, dpi=300, bbox_inches='tight')
    print("results plotted in file " + plotFileName)

