import math
import numpy as np
import os
import matplotlib.pyplot as plt
import pdb

from mpl_settings import set_plot_params


outFileName = 'jacobiSolve.dat'
vizFileName = 'jacobiSolve.jpg'
cmapChoice = 'seismic'
fwid = 10.0
fhgt = 10.0

home = os.environ['HOME']
outFilePath = os.path.join(home, 'code', 'intelShootout', 'data', outFileName)
vizFilePath = os.path.join(home, 'code', 'intelShootout', 'data', vizFileName)

with open(outFilePath, 'rb') as solutionFile:
    solutionData = solutionFile.readlines()

headerWords1 = solutionData[0].rstrip('\n').split()
headerWords2 = solutionData[1].rstrip('\n').split()
y0 = float(headerWords1[0])
yN = float(headerWords1[1])
Ny = int(headerWords1[2])
y = np.linspace(y0, yN, Ny)
x0 = float(headerWords2[0])
xN = float(headerWords2[1])
Nx = int(headerWords2[2])
x = np.linspace(x0, xN, Nx)
solution = np.zeros((Ny, Nx))
fhgt = fwid*((yN - y0)/(xN - x0))

for rowNum, row in enumerate(solutionData[2:]):
    words = row.rstrip('\n').split()
    for colNum, col in enumerate(words):
        solution[rowNum, colNum] = float(col)

vmin = abs(np.min(solution))
vmax = abs(np.max(solution))
vmax = max(vmin, vmax)
vmin = -vmax

plt.figure(1, figsize=(fwid, fhgt))
plt.pcolor(x, y, solution, vmin=vmin, vmax=vmax, cmap=cmapChoice)
plt.savefig(vizFilePath, dpi=300)
