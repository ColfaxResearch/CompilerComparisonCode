
#!/usr/bin/env python
import os
try:
    os.environ['DISPLAY']
except KeyError:
    import matplotlib as mpl
    mpl.use('Agg')
    display = False
else:
    display = True
import os
import math
import numpy as np
import argparse
from itertools import cycle
import matplotlib.pyplot as plt
from matplotlib import gridspec
import brewer2mpl
import pdb

import timingData
from mpl_settings import set_plot_params


set_plot_params(useTex=True, fontsize=16)
nProgs = 11
bmap = brewer2mpl.get_map('Set1', 'Qualitative', min(nProgs, 9))

class suitePlotter(object):
    lines = ['-','--','-.',':']

    def __init__(self, folder, logY=True, fwid=16, fhgt=10, dpi=300):
        self.linecycler = cycle(self.lines)
        self.logY = logY
        self.fwid = fwid
        self.fhgt = fhgt
        self.dpi = dpi
        self.dataFiles = list()
        if os.path.isdir(folder):
            self.folder = os.path.abspath(folder)
        else:
            raise ValueError('Invalid directory!')
            sys.exit(-1)
        self.filePaths = [os.path.join(self.folder, f) for f in os.listdir(self.folder) if os.path.isfile(os.path.join(self.folder, f))]
        self.dataFilePaths = [f for f in self.filePaths if f[-3:] == 'csv']
        for filePath in self.dataFilePaths:
            self.dataFiles.append(timingData.timingData(filePath))
        self.compilers = list(set([dataFile.compiler for dataFile in self.dataFiles]))
        self.algorithms = list(set([dataFile.algorithm for dataFile in self.dataFiles]))
        self.architectures = list(set([dataFile.architecture for dataFile in self.dataFiles]))
        self.sizes = list()
        for dataFile in self.dataFiles:
            self.sizes += dataFile.data[:, 0].tolist()
        self.sizes = list(set(self.sizes))
        self.sizes.sort()
        self.data = np.empty((len(self.architectures), len(self.algorithms), len(self.compilers), len(self.sizes), 3))
        self.data[:] = np.NAN
        self.compileTimes = np.empty(len(self.compilers))
        self.compileCounts = np.empty(len(self.compilers))
        self.suite = self.dataFiles[0].suite
        for dataFile in self.dataFiles:
            architectureNumber = self._matchArchitecture(dataFile)
            algorithmNumber = self._matchAlgorithm(dataFile)
            compilerNumber = self._matchCompiler(dataFile)
            numSizes = dataFile.data.shape[0]
            for sizeIndex in xrange(numSizes):
                sizeNumber = self._matchSize(dataFile.data[sizeIndex, 0])
                self.data[architectureNumber, algorithmNumber, compilerNumber, sizeNumber, 0] = dataFile.data[sizeIndex, 1]
                self.data[architectureNumber, algorithmNumber, compilerNumber, sizeNumber, 1] = dataFile.data[sizeIndex, 2]
                self.data[architectureNumber, algorithmNumber, compilerNumber, sizeNumber, 2] = dataFile.data[sizeIndex, 3]
            self.compileTimes[compilerNumber] += dataFile.compileTime
            self.compileCounts[compilerNumber] += 1
        self.compileTimes[:] /= self.compileCounts[:]

    def _matchAlgorithm(self, dataFile):
        for algorithmNumber, algorithm in enumerate(self.algorithms):
            if dataFile.algorithm == algorithm:
                break
        return algorithmNumber

    def _matchCompiler(self, dataFile):
        for compilerNumber, compiler in enumerate(self.compilers):
            if dataFile.compiler == compiler:
                break
        return compilerNumber

    def _matchArchitecture(self, dataFile):
        for architectureNumber, architecture in enumerate(self.architectures):
            if dataFile.architecture == architecture:
                break
        return architectureNumber

    def _matchSize(self, dataFileSize):
        for sizeNumber, size in enumerate(self.sizes):
            if dataFileSize == size:
                break
        return sizeNumber

    def plot(self):
        fig1 = plt.figure(1, figsize=(self.fwid, self.fhgt))
        numRows = 1000
        numCols = numRows
        gs = gridspec.GridSpec(numRows, numCols)
        ax1 = fig1.add_subplot(gs[:, :])
        for arch, architecture in enumerate(self.architectures):
            linestyle = next(self.linecycler)
            for a, algorithm in enumerate(self.algorithms):
                algo = algorithm.split('_')[1]
                for c, compiler in enumerate(self.compilers):
                    ax1.scatter(np.log2(self.sizes), self.data[arch, a, c,:,0], c=bmap.hex_colors[c], s=2*self.data[arch, a, c,:,2], zorder=3*c)
                    ax1.errorbar(np.log2(self.sizes), self.data[arch, a, c,:,0], color=bmap.hex_colors[c], ls=linestyle, zorder=1*c)
                    ax1.errorbar(np.log2(self.sizes), self.data[arch, a, c,:,0], yerr=self.data[arch, a, c,:,1], fmt='o', capsize=0, ls=linestyle,
                                 label='%s %s %s'%(compiler, architecture, algo), color=bmap.hex_colors[c], zorder=2*c)
        ax1.set_ylabel(r'GFLOP/s')
        ax1.set_xlabel(r'Problem Size')
        ax1.set_title(r'%s'%(self.suite))
        ax1.legend(loc=4, ncol=2, fancybox=True, fontsize=14)
        fig1.canvas.draw()
        ax1.set_xticks(np.log2(self.sizes))
        xtkLabels =  ax1.get_xticklabels()
        for index, val in enumerate(np.log2(self.sizes)):
            xtkLabels[index].set_text(r'$2^{%s}$'%(int(np.log2(self.sizes)[index])))
        ax1.set_xticklabels(xtkLabels)
        if self.logY:
            ax1.set_yscale('log')
        '''
        if self.logY:
            colStart = 675
            rowStart = 450
        else:
        '''
        colStart = 75
        rowStart = 50
        numRows = 400
        numCols = numRows
        axm1 = fig1.add_subplot(gs[rowStart:rowStart+numRows, colStart:colStart+numCols])
        width = 0.25
        axm1.set_title(r'Compile times (lower is better)', fontsize=12)
        axm1.bar(range(len(self.compileTimes)), self.compileTimes, color=bmap.hex_colors[0:len(self.compilers)])
        axm1.set_xticklabels([u''] + self.compilers, rotation=45, fontsize=10)
        axm1.set_ylabel(r'Compile time (s)', fontsize=10)
        fig1.savefig(os.path.join(self.folder, self.suite + '.jpg'), dpi=self.dpi)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process arguments for suiteRunner')
    parser.add_argument('folder', nargs=1, default=None, help='Folder to plot data from')
    parser.add_argument('--logY', dest='logY', action='store_true', help='Apply log10 to y-axis')
    parser.add_argument('--no-logY', dest='logY', action='store_false', help='Do not apply log10 to y-axis')
    parser.set_defaults(logY=True)
    parser.add_argument('--fwid', type=float, default=16.0, help='Figure width (inches)')
    parser.add_argument('--fhgt', type=float, default=10.0, help='Figure height (inches)')
    parser.add_argument('--dpi', type=float, default=300.0, help='Figure resolution (dots per inch)')
    args = parser.parse_args()
    ex = suitePlotter(folder=args.folder[0], logY=args.logY, fwid=args.fwid, fhgt=args.fhgt, dpi=args.dpi)
    ex.plot()
