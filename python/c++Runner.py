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

import math
import numpy as np
import subprocess
import shlex
import datetime
import argparse
import tempfile
import git
import psutil
import re
import warnings
import matplotlib.pyplot as plt
from matplotlib import gridspec
import brewer2mpl
import pdb

from mpl_settings import set_plot_params


set_plot_params(useTex=True, fontsize=16)
dpi = 300
fwid = 16
fhgt = 10
nProgs = 11
bmap = brewer2mpl.get_map('Set1', 'Qualitative', min(nProgs, 9))

def _parseCompilers(compilers):
    Intel = ['intel', 'icc', 'icpc', 'intel compiler', 'icc compiler', 'icpc compiler']
    GNU = ['gnu', 'gcc', 'g++', 'gnu compiler', 'gcc compiler', 'g++ compiler']
    LLVM = ['llvm', 'clang', 'clang++', 'llvm compiler', 'clang compiler', 'clang++ compiler']
    PGI = ['pgi', 'pgcc', 'pgc++', 'pgi compiler', 'pgcc compiler', 'pgc++ compiler']
    AOCC = ['amd', 'aocc', 'amd compiler', 'aocc compiler']
    ZAPCC = ['zapcc', 'zapcc++', 'zapcc compiler', 'zapcc++ compiler']
    retVal = list()
    for compiler in compilers:
        if compiler.lower() in Intel:
            retVal.append('icpc')
        elif compiler.lower() in GNU:
            retVal.append('g++')
        elif compiler.lower() in LLVM:
            retVal.append('clang++')
        elif compiler.lower() in PGI:
            retVal.append('pgc++')
        elif compiler.lower() in AOCC:
            retVal.append('aocc++')
        elif compiler.lower() in ZAPCC:
            retVal.append('zapcc++')
    return list(set(retVal))

class suiteRunner(object):
    home = os.environ['HOME']
    basedir = os.path.join(os.environ['PWD'], '..')
    res_basedir = os.path.join(os.environ['PWD'], '..', 'paper', 'results')
    cpp_basedir = os.path.join(os.environ['PWD'], '..', 'c++')
    julia_basedir = os.path.join(os.environ['PWD'], '..','julia')

    def __init__(self, suite=None, algorithm=None, architecture=None, compilers=None, sizes=None, maxthreads=None, assembly=None, runs=1, trials=10, problems=100, postfix=None):
        if assembly:
            self.assembly = assembly
        else:
            self.assembly=False
        if maxthreads:
            self.maxthreads = maxthreads
        else:
            self.maxthreads = psutil.cpu_count(logical=True)
        self.physicalcores = psutil.cpu_count(logical=False)
        self.logicalcores = psutil.cpu_count(logical=True)
        if self.maxthreads <= self.physicalcores:
            self.threadslist = np.arange(1, self.maxthreads + 1).tolist()
        elif self.maxthreads > self.physicalcores:
            self.threadslist = np.arange(1, self.physicalcores + 1).tolist()
            self.threadslist += (self.physicalcores*np.arange(2, (self.maxthreads//self.physicalcores) + 1)).tolist()
        self.repo = git.Repo(search_parent_directories=True)
        self.sha = self.repo.head.object.hexsha
        if postfix:
            if postfix == "datetime":
                self.postfix = datetime.datetime.now().isoformat()
            else:
                self.postfix = postfix
        else:
            self.postfix = self.sha[0:7]
        self.algorithm = algorithm
        self.architecture = architecture
        self.runs = runs
        self.trials = trials
        self.problems = problems
        if suite:
            self.suite = suite
        else:
            self.suite = 'sfComp'
        if compilers:
            self.compilers = _parseCompilers(compilers)
        else:
            self.compilers = ['icpc', 'g++', 'clang++', 'pgc++', 'aocc++', 'zapcc++']
        self.compileTime=len(self.compilers)*[0.0]
        if sizes:
            self.sizes = [int(size) for size in sizes]
        else:
            self.sizes = ['256', '2048']
        self.result = np.zeros((len(self.compilers), len(self.sizes), 3))

    def run(self):
        new_result = np.zeros((len(self.threadslist), 2))
        trial_result = np.zeros((self.runs, 2))
        dumpFileDirectory = os.path.join(self.res_basedir, self.suite)
        dumpFilePath = os.path.join(dumpFileDirectory, self.suite + '.' + self.postfix + '.dat')
        if not os.path.exists(dumpFileDirectory):
            os.makedirs(dumpFileDirectory)
        with open(dumpFilePath, 'ab', 0) as dump_file: # un-buffered
            for c, compiler in enumerate(self.compilers):
                clearCompileLine = 'cd %s && source %s && scons -c hbm=False auto=True architecture=%s algorithm=%s compiler=%s && cd %s || cd %s\n'%(os.path.join(self.cpp_basedir, self.suite), '~/.' + compiler + '_compiler', self.architecture, self.algorithm, compiler, self.basedir, self.basedir)
                dump_file.write(clearCompileLine)
                procClearCompile = subprocess.Popen(clearCompileLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                dump_text, err_text = procClearCompile.communicate()
                dump_file.write(dump_text)
                dump_file.write(err_text)

                if compiler == 'pgc++':
                    compileLine = 'cd %s && source %s && %s -V && time scons assembly=%s hbm=False auto=True architecture=%s algorithm=%s compiler=%s && cd %s || cd %s\n'%(os.path.join(self.cpp_basedir, self.suite), '~/.' + compiler + '_compiler', compiler, self.assembly, self.architecture, self.algorithm, compiler, self.basedir, self.basedir)
                else:
                    compileLine = 'cd %s && source %s && %s -v && time scons assembly=%s hbm=False auto=True architecture=%s algorithm=%s compiler=%s && cd %s || cd %s\n'%(os.path.join(self.cpp_basedir, self.suite), '~/.' + compiler + '_compiler', compiler, self.assembly, self.architecture, self.algorithm, compiler, self.basedir, self.basedir)
                dump_file.write(compileLine)
                procCompile = subprocess.Popen(compileLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                dump_text, err_text = procCompile.communicate()
                dump_file.write(dump_text)
                dump_file.write(err_text)

                realT = re.findall('real\t\d+m\d+.\d+s', err_text)[0]
                minT = re.findall('\d+m', realT)[0][:-1]
                secT = re.findall('\d+\.\d+s', realT)[0][:-1]
                self.compileTime[c] = float(minT)*60.0 + float(secT)

                checkLDDLine = 'cd %s && source %s && ldd %s && cd %s || cd %s\n'%(os.path.join(self.cpp_basedir, self.suite), '~/.' + compiler + '_compiler', self.suite, self.basedir, self.basedir)
                dump_file.write(checkLDDLine)
                procCheckLDD = subprocess.Popen(checkLDDLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                dump_text, err_text = procCheckLDD.communicate()
                dump_file.write(dump_text)
                dump_file.write(err_text)

                for s, size in enumerate(self.sizes):
                    notTurnedOver = True
                    exceededCores = False
                    for t, threads in enumerate(self.threadslist):
                        for n in xrange(self.runs):
                            if compiler == 'pgc++':
                                runSuiteLine = 'cd %s && source %s && OMP_PROC_BIND=true OMP_NUM_THREADS=%s PROBLEM_SIZE=%s NUM_PROBLEMS=%s NUM_TRIALS=%s ./%s\n'%(os.path.join(self.cpp_basedir, self.suite), '~/.' + compiler + '_compiler', threads, size, self.problems, self.trials, self.suite)
                            else:
                                runSuiteLine = 'cd %s && source %s && OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=%s PROBLEM_SIZE=%s NUM_PROBLEMS=%s NUM_TRIALS=%s ./%s\n'%(os.path.join(self.cpp_basedir, self.suite), '~/.' + compiler + '_compiler', threads, size, self.problems, self.trials, self.suite)
                            dump_file.write(runSuiteLine)
                            procRunSuite = subprocess.Popen(runSuiteLine, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                            run_text, err_text = procRunSuite.communicate()
                            dump_file.write(run_text)
                            dump_file.write(err_text)
                            trial_result[n,0] = float(run_text.rstrip('\n').split()[2])
                            trial_result[n,1] = float(run_text.rstrip('\n').split()[4])
                        new_result[t, 0] = np.max(trial_result[:,0])
                        new_result[t, 1] = trial_result[np.where(new_result[t,0] == trial_result[:,0])[0][0], 1]
                    self.result[c,s,0] = np.max(new_result[:, 0])
                    bestthreads = np.where(new_result[:,0] == self.result[c,s,0])[0][0]
                    self.result[c,s,1] = new_result[bestthreads, 1]
                    self.result[c,s,2] = self.threadslist[bestthreads]
                    new_result[:, 0] = 0.0
                    new_result[:, 1] = 0.0

    def write(self):
        for s, size in enumerate(self.sizes):
            for c, compiler in enumerate(self.compilers):
                write_file_path = os.path.join(self.res_basedir, self.suite, self.suite + '.' + compiler + '.%s.%s.%s.%s.'%(size, self.runs, self.trials, self.problems) + self.postfix + '.csv')
                with open(write_file_path, 'wb', 0) as file_path: # un-buffered
                    header = '# compile time = %17.16e s\n'%(self.compileTime[c])
                    file_path.write(header)
                    header = '# %12s %22s %22s %12s\n'%('PROBLEM_SIZE', 'GFLOPS/s','ERR_GFLOPS/s','NUM_THREADS')
                    file_path.write(header)
                    line = '%14d %17.16e %17.16e %12d\n'%(size, self.result[c, s, 0], self.result[c, s, 1], self.result[c, s, 2])
                    file_path.write(line)

    def plot(self):
        vmin = abs(np.min(self.result))
        vmax = abs(np.max(self.result))
        vmax = max(vmin, vmax)
        vmin = -vmax
        fig1 = plt.figure(1, figsize=(fwid, fhgt))
        numRows = 1000
        numCols = numRows
        gs = gridspec.GridSpec(numRows, numCols)
        ax1 = fig1.add_subplot(gs[:, :])
        for c, compiler in enumerate(self.compilers):
            ax1.scatter(np.log2(self.sizes), self.result[c,:,0], c=bmap.hex_colors[c], s=self.result[c,:,2]**2)
            ax1.errorbar(np.log2(self.sizes), self.result[c,:,0], color=bmap.hex_colors[c])
            ax1.errorbar(np.log2(self.sizes), self.result[c,:,0], yerr=self.result[c,:,1], fmt='o', capsize=0,
                         label=compiler, color=bmap.hex_colors[c])
        ax1.set_xlabel(r'Problem Size')
        ax1.set_ylabel(r'GFLOP/s')
        ax1.set_title(r'Structure Function computation')
        ax1.legend(loc=4, ncol=2, fancybox=True, fontsize=14)
        fig1.canvas.draw()
        ax1.set_xticks(np.log2(self.sizes))
        tkLabels =  ax1.get_xticklabels()
        for index, val in enumerate(np.log2(self.sizes)):
            tkLabels[index].set_text(r'$2^{%s}$'%(int(np.log2(self.sizes)[index])))
        ax1.set_xticklabels(tkLabels)
        fig1.savefig(os.path.join(self.res_basedir, self.suite + '.' + self.postfix + '.jpg'), dpi=dpi)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process arguments for suiteRunner')
    parser.add_argument('suite', nargs=1, default=None, help='Suite to run')
    parser.add_argument('algorithm', nargs=1, default=None, help='Algorithm to run')
    parser.add_argument('architecture', nargs=1, default=None, help='Architecture to run on')
    parser.add_argument('--compilers', nargs='*', default=None, help='List of compilers to run')
    parser.add_argument('--sizes', nargs='*', default=None, help='List of problem sizes to run')
    parser.add_argument('--postfix', default=None, help='Postfix to be applied to output files')
    parser.add_argument('--maxthreads', type=int, default=None, help='Maximum number of threads to test')
    parser.add_argument('--assembly', type=bool, default=False, help='Output assembly?')
    parser.add_argument('-r', '--runs', type=int, default=1, help='Number of times to run each suite')
    parser.add_argument('-t', '--trials', type=int, default=10, help='Number of trials per suite')
    parser.add_argument('-p', '--problems', type=int, default=100, help='Number of problems (vectors/matrices etc...) per suite')
    args = parser.parse_args()
    if args.trials < 3:
        warnings.warn('The first 3 trials are for warm-up & are usually discarded! Expect the performance to be 0.0 GFLOPS/s')

    ex = suiteRunner(suite=args.suite[0], algorithm=args.algorithm[0], architecture=args.architecture[0], compilers=args.compilers, sizes=args.sizes, maxthreads=args.maxthreads, assembly=args.assembly, runs=args.runs, trials=args.trials, problems=args.problems, postfix=args.postfix)
    ex.run()
    ex.write()
    #if display:
        #ex.plot()
