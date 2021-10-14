import cmd #for program commands
import argparse #for figure arguments
import numpy as np #for... need I explain?
import matplotlib #urm... seems kinda obvious
import matplotlib.pyplot as plt #convenience
import sys
import os

from .data import xyData, stickData, cmgData, duoOutData, roueffData, detect_filetype
from .view import GraphView, StickView, linelistComparisonView

### User interface object

class cmgPrompt(cmd.Cmd):
    """The cmdGraph user command prompt interface object.

    The cmdPrompt controls the operation of the program for a given figure.
    In normal operation, with a single figure instance, this amounts to the 
    entire program. 
    
    The class defines the program commands, reads/writes data files, and owns
    and passes commands to the current View instance.
    """
    prompt = '> ' #set left hand icon for prompt

    def __init__(self, mode='launch', **kwargs):
        """Initialise a new prompt mode. 
        
        Creates a figure and sets the View object and Data type specified by the 
        'mode' keyword argument. The default value is 'launch', which is used 
        when the program starts to call the Cmd superclass __init__ function. 
        Other keyword arguments are passed to the Cmd superclass.
        """
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.ion()
        self._modes = ['graph', 'stick', 'linelistComparison'] #define available modes
        if mode == 'launch': #for first __init__
            cmd.Cmd.__init__(self, **kwargs)
            mode = 'stick' #default to stick
        if mode in self._modes:
            self.mode = mode
            self.fig = plt.figure()
            self._single = False
            if mode == 'graph':
                self._View = GraphView(self.fig) #set the view mode
            elif mode == 'stick':
                self._View = StickView(self.fig)
            elif mode == 'linelistComparison':
                self._View = linelistComparisonView(self.fig)
        else:
            print("Mode not defined.")

    def onecmd(self, line):
        """Overwritten from default cmd.onecmd method to prevent exiting command
        loop when an error is raised.

        """
        cmd, arg, line = self.parseline(line)
        if not line:
            return self.emptyline()
        if cmd is None:
            return self.default(line)
        self.lastcmd = line
        if line == 'EOF' :
            self.lastcmd = ''
        if cmd == '':
            return self.default(line)
        else:
            try:
                func = getattr(self, 'do_' + cmd)
                return func(arg)
            except AttributeError:
                return self.default(line)

    ## Command methods ##
    def default(self, line):
        """Pass input string to the View instance's argparse parser if the string
        is unrecognised as a cmdPrompt command. Check if argparse has successfully
        parsed and warn if not.
        
        """
        argIn = line.split(' ') #argparse expects list of arguments
        hasParsed = self._View.parse(argIn)
        if not hasParsed: #unrecognised argparse argument
            print("Unrecognised argument.")
        else:
            self._changed = True

    def do_exit(self, inp): #not working
        """Checks if in single plot mode, if so makes all View plots live,
        and exits single mode. Otherwise, asks before exiting program.
        
        """
        if self._single:
            self._View.livePlots = self._View.plots
            self._single = False
            print("Leaving single plot mode...")
            return False
        else:
            yn = 'm'
            while yn not in ['Y', 'y', 'N', 'n']:
                yn = input("Exit the program? [y/n] ")
            if yn in ['N', 'n']:
                print("Exit cancelled.")
                return False
            elif yn in ['Y', 'y']:
                print("Exiting cmdGraph...")
                return True          
    def help_exit(self):
        print("usage: exit\n    Exit current environment, or if no environment, exit program.")

    def do_adat(self, inp):
        """Calls the read file method to initialise a Data instance for the input
        file using the active data mode. Requires that self._Data has been set 
        to a valid Data class by __init__. Then adds Data to the View instance as
        as plot.
        
        """
        for inFile in inp.split(): #assume spaces in input delimit files
            fileType = detect_filetype(inFile) #detect data type
            _Data = fileType(inFile).read_file() #read input file and return data object
            self._View.add_plot(_Data) #call view mode add plot function
    def help_adat(self):
        print("usage: adat <file1> <file2> ... \n    Add lines(s) to figure from file(s). ")

    def do_ddat(self, inp):
        """Remove a data file from the figure. Not working yet, placeholder only.
        
        """
        for plot in self._View.plots:
            if plot._Data.fName == inp:
                self._View.plots.remove(plot)
    def help_ddat(self):
        print("usage: ddat <file>\n    Placeholder.")

    def do_single(self, inp):
        """Replace live plots in View instance with a single plot, allows for
        tweaking of individual Plot instances within a figure without having to
        pass arguments for every other instance.
        
        Requires self.do_exit() call to terminate and return all plots to View
        live plots.
        
        """
        for plot in self._View.plots:
            _Data = plot._Data
            if _Data.fName == inp:
                self._View.livePlots = [plot] #single live plot
                self._single = True #enter single mode
                print("Entering single plot mode for '{}'".format(_Data.fName))
                break
        if not self._single:
            print("Plot file not found.")
    def help_single(self):
        print("usage: single <file>\n    Switch to single dataset edit mode.")

    def do_mode(self, inp):
        """Set the current plotting mode. Calls the cmdPrompt __init__ method
        to initialise a new View instance and change the data mode. Closes current
        figure (without saving) and opens new figure for the input View mode.
        
        """
        if inp in self._modes:
            plt.close(fig='all')
            self.__init__(mode=inp)
        else:
            print("Mode not defined.")
    def help_mode(self):
        print("usage: mode <mode>\n    Change figure mode (graph, stick)")

    def do_save(self, inp):
        """Save configuration for current figure in a native cmdGraph 'save file'.
        This is done by initialising an instance o the cmdGraph data class and
        calling it's write file method for the current View instance.
        
        """
        if not inp:
            if not hasattr(self, '_savefile'):
                i = 1
                while os.path.isfile("./autosave{0}.cmg".format(i)):
                    i +=1
                self._savefile = "autosave{0}.cmg".format(i)
            inp = self._savefile
        cmgData(inp).write_file(self._View)
        self._changed = False
        self._savefile = inp
    def help_save(self):
        print("usage: save <filename>\n    Save the current figure as <filename>.cmg, allows figures to be reload after exit.")

    def do_load(self, inp):
        """Load a View instance from a cmdGraph 'save file' by initialising a
        cmdGraph data instance and calling it's read file method. Then pass
        each item in the data attribute as though it were a user command.
        """
        buff = cmgData(inp).read_file()
        for cmdString in buff.dat:
            self.onecmd(cmdString)
        if self._single:
            self.do_exit('')
        self._savefile = inp
    def help_load(self):
        print("usage: load <filename>.cmg\n    Load figure from cmdGraph file.")

    def do_print(self, inp):
        """Rudimentary print figure to file. Redirects to plt.savefig(), the
        filetpye is automatically chosen from the suffix used in the input
        string, i.e '.pdf', '.jpg'.
        """
        if not inp:
            if not hasattr(self, '_printfile'):
                i = 1
                while os.path.isfile("./autoprint{0}.pdf".format(i)):
                    i +=1
                self._printfile = "autoprint{0}.pdf".format(i)
            inp = self._printfile
        plt.savefig(inp)
    def help_print(self):
        print("usage print <filename>.<filetype>\n    Print figure to file.")

    def do_tight(self, inp):
        """Apply matplotlib tight_layout."""
        plt.tight_layout()
    def help_tight(self):
        print("usage: tight\n    Apply tight layout.")

def run():
    """Run the cmdGraph program."""
    cmgPrompt(mode='launch').cmdloop() #run program
