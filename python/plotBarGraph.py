import math
import numpy as np
import os
from itertools import cycle
import brewer2mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pdb

from mpl_settings import set_plot_params


fwid = 20
fhgt = 7
dpi = 300
set_plot_params(useTex=True, fontsize=26)
nProgs = 11
bmap = brewer2mpl.get_map('Pastel1', 'Qualitative', min(nProgs, 6))
hexColors = [bmap.hex_colors[2], bmap.hex_colors[0], bmap.hex_colors[5], bmap.hex_colors[3],
             bmap.hex_colors[4], bmap.hex_colors[1]]

pythonDir = os.environ['PWD']
resultsDir = os.path.join(pythonDir, '..', 'paper', 'results')

arch = 'bdw'

sfCompOBLKIO_VEC = 'sfComp_' + arch + '_OBLKIO_VEC'
sfCompOBLKIO_OPT = 'sfComp_' + arch + '_OBLKIO_OPT'
luDecompKIJ_VEC = 'luDecomp_' + arch + '_KIJ_VEC'
luDecompKIJ_OPT = 'luDecomp_' + arch + '_KIJ_OPT'
jacobiSolveSOLVE_VEC = 'jacobiSolve_' + arch + '_SOLVE_VEC'
jacobiSolveSOLVE_OPT = 'jacobiSolve_' + arch + '_SOLVE_OPT'
tmvCompile = 'tmv_COMPILE'
sfCompOBLKIO_VECPath = os.path.join(resultsDir, sfCompOBLKIO_VEC)
sfCompOBLKIO_OPTPath = os.path.join(resultsDir, sfCompOBLKIO_OPT)
luDecompKIJ_VECPath = os.path.join(resultsDir, luDecompKIJ_VEC)
luDecompKIJ_OPTPath = os.path.join(resultsDir, luDecompKIJ_OPT)
jacobiSolveSOLVE_VECPath = os.path.join(resultsDir, jacobiSolveSOLVE_VEC)
jacobiSolveSOLVE_OPTPath = os.path.join(resultsDir, jacobiSolveSOLVE_OPT)
tmvCompilePath = os.path.join(resultsDir, 'tmv_COMPILE')
vecResultNames = [luDecompKIJ_VEC, jacobiSolveSOLVE_VEC, sfCompOBLKIO_VEC]
vecResultPaths = [luDecompKIJ_VECPath, jacobiSolveSOLVE_VECPath, sfCompOBLKIO_VECPath]
optResultNames = [luDecompKIJ_OPT, jacobiSolveSOLVE_OPT, sfCompOBLKIO_OPT]
optResultPaths = [luDecompKIJ_OPTPath, jacobiSolveSOLVE_OPTPath, sfCompOBLKIO_OPTPath]
timeResultNames = [tmvCompile]
timeResultPath = [tmvCompilePath]
compilerList = ['PGC++', 'Clang', 'AOCC', 'Zapcc', 'G++', 'Intel C++']
hatchPatterns = ['+', '\\', '-', 'x', '|', '/']
hatchCycler = cycle(hatchPatterns)
lineStyles = ['solid', 'dashed', 'dashdot', 'dotted', '-', '--', '-.', ':', 'None', ' ', '']
lineCycler = cycle(lineStyles)


def vecBarPlot(resultArray):
    fig = plt.figure(1, figsize=(fwid, fhgt))
    numRows = 1000
    numCols = numRows
    gs = gridspec.GridSpec(numRows, numCols)
    ax1 = fig.add_subplot(gs[:, :])
    barWidth = 0.05
    ax1.axhline(1.0, color='#989898', linewidth=1, linestyle='--', alpha = 25, zorder=0)
    for res in range(resultArray.shape[0]):
        for c, compiler in enumerate(compilerList):
            speedUp = resultArray[res, c, 1]/resultArray[res, 4, 1]
            if res == 0:
                ax1.bar(9*res*barWidth + (c - 3)*barWidth, speedUp,
                        barWidth, align='edge', color=hexColors[c], hatch=hatchPatterns[c],
                        edgecolor='#000000', linewidth=1.0, ls='solid', label=compiler, zorder=5)
            else:
                ax1.bar(9*res*barWidth + (c - 3)*barWidth, speedUp,
                        barWidth, align='edge', color=hexColors[c], hatch=hatchPatterns[c],
                        edgecolor='#000000', linewidth=1.0, ls='solid', zorder=5)
            ax1.text(9*res*barWidth + (c - 3)*barWidth + barWidth/2, speedUp + 0.01, r'$%3.2f$'%(speedUp),
                    ha='center', va='bottom', fontsize=20)
    ax1.set_xticks([0.0, 9.0*barWidth, 18.0*barWidth])
    ax1.set_xticklabels([r'LU Decomposition', r'Jacobi Solver', r'SF Computation'])
    ax1.set_ylabel(r'Relative Performance (G++ $ = 1$)')
    plt.legend(bbox_to_anchor=(0.025, 0.85, 0.95, .102), loc=3, ncol=6, mode="expand", borderaxespad=0.)
    plt.ylim(ymax=2.0)
    fig.savefig(os.path.join(resultsDir, arch + '_VEC.pdf'))
    fig.savefig(os.path.join(resultsDir, arch + '_VEC.jpg'))
    fig.clf()


def optBarPlot(resultArray):
    fig = plt.figure(2, figsize=(fwid, fhgt))
    numRows = 1000
    numCols = numRows
    gs = gridspec.GridSpec(numRows, numCols)
    ax1 = fig.add_subplot(gs[:, :])
    barWidth = 0.05
    ax1.axhline(1.0, color='#989898', linewidth=1, linestyle='--', alpha = 25, zorder=0)
    for res in range(resultArray.shape[0]):
        for c, compiler in enumerate(compilerList):
            speedUp = resultArray[res, c, 1]/resultArray[res, 4, 1]
            if res == 0:
                ax1.bar(9*res*barWidth + (c - 3)*barWidth, speedUp,
                        barWidth, align='edge', color=hexColors[c], hatch=hatchPatterns[c],
                        edgecolor='#000000', linewidth=1.0, ls='solid', label=compiler)
                ax1.text(9*res*barWidth + (c - 3)*barWidth + barWidth/2, speedUp + 0.01, r'$%3.2f$'%(speedUp),
                         ha='center', va='bottom', fontsize=20)
            else:
                if res == 2 and c == 0:
                    ax1.bar(9*res*barWidth + (c - 3)*barWidth, np.nan,
                            barWidth, align='edge', color=hexColors[c], hatch=hatchPatterns[c],
                            edgecolor='#000000', linewidth=1.0, ls='solid')
                else:
                    ax1.bar(9*res*barWidth + (c - 3)*barWidth, speedUp,
                            barWidth, align='edge', color=hexColors[c], hatch=hatchPatterns[c],
                            edgecolor='#000000', linewidth=1.0, ls='solid')
                    ax1.text(9*res*barWidth + (c - 3)*barWidth + barWidth/2, speedUp + 0.01, r'$%3.2f$'%(speedUp),
                             ha='center', va='bottom', fontsize=20)
    ax1.set_xticks([0.0, 9.0*barWidth, 18.0*barWidth])
    ax1.set_xticklabels([r'LU Decomposition', r'Jacobi Solver', r'SF Computation'])
    ax1.set_ylabel(r'Relative Performance (G++ $ = 1$)')
    plt.legend(bbox_to_anchor=(0.025, 0.85, 0.95, .102), loc=3, ncol=6, mode="expand", borderaxespad=0.)
    plt.ylim(ymax=2.25)
    fig.savefig(os.path.join(resultsDir, arch + '_OPT.pdf'))
    fig.savefig(os.path.join(resultsDir, arch + '_OPT.jpg'))
    fig.clf()


def timeBarPlot(resultArray):
    fig = plt.figure(3, figsize=(fwid/2.1, fwid/2.1))
    numRows = 1000
    numCols = numRows
    gs = gridspec.GridSpec(numRows, numCols)
    ax1 = fig.add_subplot(gs[:, :])
    barWidth = 0.05
    ax1.axhline(1.0, color='#989898', linewidth=1, linestyle='--', alpha = 25, zorder=0)
    for res in range(resultArray.shape[0]):
        for c, compiler in enumerate(compilerList):
            speedUp = resultArray[res, 4, 1]/resultArray[res, c, 1]
            ax1.bar(9*res*barWidth + (c - 3)*barWidth, speedUp,
                    barWidth, align='edge', color=hexColors[c], hatch=hatchPatterns[c],
                    edgecolor='#000000', linewidth=1.0, ls='solid', label=compiler)
            ax1.text(9*res*barWidth + (c - 3)*barWidth + barWidth/2, speedUp + 0.01, r'$%3.2f$'%(speedUp),
                    ha='center', va='bottom', fontsize=20)
    ax1.set_xticks([0.0])
    ax1.set_xticklabels([r'tmv Compile'])
    ax1.set_ylabel(r'Relative Compile Speed (G++ $ = 1$)')
    plt.legend(bbox_to_anchor=(0.025, 0.825, 0.95, .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
    plt.ylim(ymax=2.25)
    fig.savefig(os.path.join(resultsDir, 'TIME.pdf'))
    fig.savefig(os.path.join(resultsDir, 'TIME.jpg'))
    fig.clf()


resultArray = np.zeros((3, len(compilerList), 5))

for res, resultPath in enumerate(vecResultPaths):
    with open(resultPath + '.dat', 'rb') as resultFile:
        resultLines = resultFile.readlines()
    for compiler, resultLine in enumerate(resultLines):
        resultWords = resultLine.rstrip('\n').split()
        resultArray[res, compiler, 0] = int(resultWords[1])
        resultArray[res, compiler, 1] = float(resultWords[2])
        resultArray[res, compiler, 2] = float(resultWords[3])
        resultArray[res, compiler, 3] = int(resultWords[4])
        resultArray[res, compiler, 4] = float(resultWords[5])
vecBarPlot(resultArray)

for res, resultPath in enumerate(optResultPaths):
    with open(resultPath + '.dat', 'rb') as resultFile:
        resultLines = resultFile.readlines()
    for compiler, resultLine in enumerate(resultLines):
        resultWords = resultLine.rstrip('\n').split()
        resultArray[res, compiler, 0] = int(resultWords[1])
        resultArray[res, compiler, 1] = float(resultWords[2])
        resultArray[res, compiler, 2] = float(resultWords[3])
        resultArray[res, compiler, 3] = int(resultWords[4])
        resultArray[res, compiler, 4] = float(resultWords[5])
optBarPlot(resultArray)

resultArray = np.zeros((1, len(compilerList), 5))

for res, resultPath in enumerate(timeResultPath):
    with open(resultPath + '.dat', 'rb') as resultFile:
        resultLines = resultFile.readlines()
    for compiler, resultLine in enumerate(resultLines):
        resultWords = resultLine.rstrip('\n').split()
        resultArray[res, compiler, 0] = int(resultWords[1])
        resultArray[res, compiler, 1] = float(resultWords[2])
        resultArray[res, compiler, 2] = float(resultWords[3])
        resultArray[res, compiler, 3] = int(resultWords[4])
        resultArray[res, compiler, 4] = float(resultWords[5])
timeBarPlot(resultArray)
