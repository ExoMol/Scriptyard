import cmd #for program commands
import argparse #for figure arguments
import numpy as np #for... need I explain?
import matplotlib #urm... seems kinda obvious
import matplotlib.pyplot as plt #convenience
import sys
import pandas as pd

_small_ = 1e-300

### Data Classes 

class _Data:
    def __init__(self, filename):
        self.fName = filename
        self.dat = None

class xyData(_Data):
    """Data format for two column x, y data files where each row is a data point.
    assumes no column headers. Currently just a skeleton class for basic operation."""
    def read_file(self):
        """Read x, y columns from file"""
        with open(self.fName, 'r') as f:
            datLines = []
            for l in f:
                line = l.strip().split()
                datLines.append(line)
            self.dat = np.array(datLines, dtype=float).T
        return self

class stickData(_Data):
    """Data format for stick spectra. Input file with two columns of x positions and
    stick height is converted into a set of lines for matplotlib LineCollection."""
    def read_file(self):
        """Read wavenumber, intensity columns from file"""
        with open(self.fName, 'r') as f:
            datLines = []
            ymin = ymax = 0
            for l in f:
                x, y = l.strip().split()
                datLines.append([[x, _small_], [x, y]])
                ymin = float(y) if (float(y) < ymin) else ymin
                ymax = float(y) if (float(y) > ymax) else ymax
            self.dat = np.array(datLines, dtype=float)
            self.ybounds = [float(ymin), (ymax)]
        return self

class cmgData(_Data):
    """Data format for internal cmdGraph data files. This allows configurations
    to be saved so that they can be transferred and reloaded after the program
    exits.

    The format of cmdGraph save files is simply a list of commands to be executed 
    to reproduce the configuration from the input data files. A standard save file
    would look like:

    | ---cmdGraph---                    <- file header confirms cmdGraph format
    | mode graph                        <- select view mode
    | adat file1.csv file2.csv          <- load data files
    | --xrange 0 10                     <- axes arguments
    | --yrange -1 1                     <-
    | single file1.csv                  <- enter single mode for first file
    | --linewidth 2.0                   <- plot arguments of first file
    | --label The#first#data
    | single file2.csv                  <- enter single mode for second file
    | --label The#second#data           <- plot arguments of second file
    | --linecolour red

    """
    def write_file(self, view):
        """Write the View instance to a '.cmg' file with provided name. 
        
        Tries to retrieve corresponding '_prop_' attribute for every argument in
        the View instance's list of argument names. First at the level of View
        attributes for axes arguments, or else at the level of Plot attributes
        for plot arguments.

        """

        f = open(self.fName, 'w')
        f.write("---cmdGraph---" + '\n')
        f.write("mode " + view.mode + '\n')
        f.write(("adat " + " {}"*len(view.plots)
            + '\n').format(*[P._Data.fName for P in view.plots]))
        _plotArgs = []
        for argName in view.argNames:
            try:
                argVal = getattr(view, '_prop_' + argName)
                argVal = str(argVal) if (type(argVal) is not str) else argVal
                f.write('--' + argName + ' ' + argVal + '\n')
            except:
                _plotArgs.append(argName)
        for plot in view.plots:
            f.write('single ' + plot._Data.fName + '\n')
            for argName in _plotArgs:
                try:
                    argVal = getattr(plot, '_prop_' + argName)
                    argVal = str(argVal) if (type(argVal) is not str) else argVal
                    f.write('--' + argName + ' ' + argVal + '\n')
                except:
                    pass
        f.close()

    def read_file(self):
        """Read a View instance from file and return the list of command strings 
        to be executed in order to recover the figure.

        """
        with open(self.fName, 'r') as f:
            datLines = []
            start = False
            for line in f:
                line = line.strip()
                while not start:
                    if line == '---cmdGraph---':
                        start = True
                        line = ''
                    continue
                if line == '':
                    continue
                else:
                    datLines.append(line)
        if not start:
            print("Not recognised as a cmdGraph file.")
        else:
            self.dat = datLines
        return self

class duoOutData(_Data):
    """Data format for reading transition data from Duo '.out' files. Automatically
    locates and extracts the transitions data from the Einstein coefficients and
    linestrengths section of the output file. Returns a numpy array of relevant values."""
    def read_file(self):
        """Read linelist from Duo output file."""
        header = "    J Gamma <-   J  Gamma Typ       Ei     <-      Ef          nu_if        S(f<-i)          A(if)            I(f<-i)       State v lambda sigma  omega <- State v lambda sigma  omega "
        with open(self.fName, 'r') as f:
            datLines = []
            atBlock = False
            l = 0
            while True:
                l += 1
                line = f.readline().rstrip('\n')
                if line == header:
                    atBlock = True
                    continue
                if atBlock and line:
                    if line == "done":
                        break
                    line = line.replace('<-', '').replace('(', '').replace(')', '')
                    datLines.append(
                        line.replace('<-', '').replace('(', '').replace(')', '').split()
                        ) #remove non-data columns
        self.dat = np.array(datLines)
        self.cols = [
            'rotational_final',
            'gamma_final',
            'rotational_initial',
            'gamma_initial',
            'transition_branch',
            'energy_final_cm',
            'energy_initial_cm',
            'wavenumber',
            'linestrength_S',
            'einstein_A',
            'intensity_I',
            'electronic_final',
            'vibrational_final',
            'lambda_final',
            'sigma_final',
            'omega_final',
            'electronic_initial',
            'vibrational_initial',
            'lambda_initial',
            'sigma_initial',
            'omega_initial'
            ]
        self.typedict = {key: float for key in self.cols}
        self.typedict.update({'gamma_initial': str, 'gamma_final': str, 'transition_branch': str})
        return self

class roueffData(_Data):
    """Data format for linelist in the format given by Roueff et al. 2019."""
    def read_file(self):
        """Read linelist from Komasa format."""
        with open(self.fName, 'r') as f:
            datLines = []
            while True:
                line = f.readline().rstrip('\n')
                if line:
                    datLines.append(line.split())
                else:
                    break
        i = input("Select E2 (0), M1 (1) or both (2) lines: ")
        choose = [8, 9, 10]
        self.dat = np.array(datLines)[:, [*range(8), choose[int(i)], *range(12,16)]]
        self.cols = [
            'vibrational_final',
            'rotational_final',
            'vibrational_initial',
            'rotational_initial',
            'wavenumber',
            'wavenumber_error',
            'wavelength',
            'wavelength_error',
            'einstein_A',
            'energy_final_cm',
            'energy_final_error',
            'energy_final_kelvin',
            'g_factor_final'
            ]
        self.typedict = {
            key: float for key in self.cols
            }
        return self
        
def detect_filetype(fName):
    pref, suf = fName.split('.')
    if suf == 'stick':
        return stickData
    elif suf == 'out':
        return duoOutData
    elif suf == 'cmg':
        return cmgData
    elif suf == 'txt':
        f = open(fName, 'r')
        firstLine = f.readline()
        if len(firstLine.split()) == 2:
            return xyData
        elif len(firstLine.split()) == 16:
            return roueffData
    elif suf == 'csv':
        return xyData
    else:
        print("File type not recognised")

