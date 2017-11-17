import math
import numpy as np
import os
import pdb


class timingData(object):
    def __init__(self, dataFilePath):
        self.headerLineCount = 0
        self.headerLines = list()
        if not os.path.isfile(dataFilePath) or dataFilePath[-3:] != 'csv':
            raise ValueError('Invalid dataFilePath: %s'%(dataFilePath))
        self.dataFilePath = dataFilePath
        folder, fileName = os.path.split(self.dataFilePath)
        parts = fileName.split('.')
        self.suite = parts[0]
        self.compiler = parts[1]
        self.runs = int(parts[3])
        self.trials = int(parts[4])
        self.problems = int(parts[5])
        self.algorithm = parts[6]
        self.architecture = parts[7]
        with open(dataFilePath, 'rb') as dataFile:
            lines = dataFile.readlines()
        for line in lines:
            words = line.rstrip('\n').split()
            if words[0] == '#':
                self.headerLineCount += 1
                self.headerLines.append(line)
                if words[1] == 'compile' and words[2] == 'time':
                    self.compileTime = float(words[4])
        self.dataLineCount = 0
        self.data = np.zeros((len(lines) - self.headerLineCount, 4))
        for line in lines:
            words = line.rstrip('\n').split()
            if words[0] != '#':
                self.data[self.dataLineCount, 0] = int(words[0])
                self.data[self.dataLineCount, 1] = float(words[1])
                self.data[self.dataLineCount, 2] = float(words[2])
                self.data[self.dataLineCount, 3] = int(words[3])
                self.dataLineCount += 1
