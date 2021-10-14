import cmd #for program commands
import argparse #for figure arguments
import numpy as np #for... need I explain?
import matplotlib #urm... seems kinda obvious
import matplotlib.pyplot as plt #convenience
import pandas as pd
import sys

from .data import xyData, stickData, cmgData, duoOutData, roueffData, detect_filetype


### View Classes

class _View:
    """The superclass for View modes. 
    
    Defines basic properties required for a View object to function within the
    program, but should always be expanded with custom methods.

    Every View instance has a Plot class. Instances of a Plot class are
    associated with a data file for the current View instance, and have a set
    of commands that can be passed by the user to edit their properties and call
    methods. View instances also have an axes.
    
    Currently axes are not defined as an inner class, but this might change if
    View modes are introduced with multiple axes.
    
    See existing subclasses for examples.
    """
    def __init__(self, figure=None, mode=None):
        """View mode superclass with necessary functions, do not overwrite."""
        matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=[
            "#2271b2",
            "#e69f00",
            "#359b73",
            "#f748a5",
            "#f0e442",
            "#3db7e9",
            "#d55e00"
        ])
        self.fig       = figure
        self.mode      = mode
        self.plots     = []
        self.livePlots = self.plots
        self.argNames  = []
        self.parser    = argparse.ArgumentParser(
            usage="Figure commands should have a '-' or '--' in front.\n"
            "       Type '-h' for a list of figure commands available in the current mode."
        )
        self._add_arg('-xl', '--xlabel', nargs='+', type=str, metavar='str')
        self._add_arg('-yl', '--ylabel', nargs='+', type=str, metavar='str')
        self._add_arg('-fs', '--figsize', nargs=2, type=str, metavar='float')
        self._add_arg('-ys', '--yscale', nargs=1, type=str, metavar='str',
                help="Set the scale of the y axis, e.g 'log', 'linear'")
        self._add_arg('-xr', '--xrange', nargs=2, type=str, metavar='float',
            help="min and max range of figure x-axis, use * * to automatically choose value")
        self._add_arg('-yr', '--yrange', nargs=2, type=str, metavar='float',
            help="min and max range of figure y-axis, use * * to automatically choose value")


    def parse(self, inp):
        """Method for parsing arguments using the parser defined for the future
        View subclass. Tries argument and if not parsed returns false to the
        cmdPrompt instance. Otherwise attempts to set attribute of axes, and
        then of plots if command is not recognised as an axes property.

        """
        try: #catch argparse system exit
            args = self.parser.parse_args(inp) #use argparse to parse args
        except:
            return False
        for argName in self.argNames:
            if hasattr(args, argName):
                inp = getattr(args, argName)
                if hasattr(self, '_set_'+argName): #axes command?
                    _set = getattr(self, '_set_'+argName)
                    _set(inp)
                else: #plot command?
                    for p, plot in enumerate(self.livePlots):
                        _set = getattr(plot, '_set_'+argName)
                        _set(inp[p])
        return True

    def _add_arg(self, *args, **kwargs):
        """Adds argument to View instance, essentially a wrapper for argparse
        'add_argument' that forces default to SUPPRESS and adds the argument
        name to an internal list.

        Every argument should have an associated 'self._set_<argument>' method,
        the first action of which should always be to set 'self._prop_<argument>'
        equal to the console input. Without '_prop_' attributes, the
        configuration cannot be saved, and without '_set_' methods the argument
        will not activate any change.
        """
        self.parser.add_argument(*args, **kwargs, default=argparse.SUPPRESS)
        self.argNames.append(args[1].lstrip('--'))

    class Plot:
        def __init__(self, data, ax):
            """Template class for the View subclass's Plot class. The Plot class
            defines the type of graph the user is plotting, and has associated
            properties and methods. All plot instances should however have an
            associated data object (data) and axes (ax).
            
            The __init__ function is called by the cmdPrompt adat method when
            a data file is added to the current View instance.
            """
            self._Data = data
            self.ax  = ax

    def add_plot(self, data):
        """Method for adding a Plot instance to the current View with data from
        a given Data object. Adds plot to list of plots in current View and resets
        live plots. Assumes plot is being added to current matplotlib axes.
        """
        plot = self.Plot(data, plt.gca()) #plot object init method creates matplotlib artist
        self.plots.append(plot) #add to list of plots
        self.livePlots = self.plots #update live plots

    def _set_xlabel(self, inp):
        """Set x label for axis."""
        self._prop_xlabel = inp[0]
        self.ax.set_xlabel(inp[0].replace('#', ' '))
    def _set_ylabel(self, inp):
        """Set y label for axis."""
        self._prop_ylabel = inp[0]
        self.ax.set_ylabel(inp[0].replace('#', ' '))
    def _set_figsize(self, inp):
        """Set figure size."""
        self._prop_figsize = "{} {}".format(*inp)
        self.fig.set_size_inches(float(inp[0]), float(inp[1]))
    def _set_yscale(self, inp):
        """Set y scale."""
        self._prop_yscale = inp[0]
        self.ax.set_yscale(inp[0])
    # Axes methods
    def _set_xrange(self, inp):
        """Set the x range of axis. If an asteriks is detected in place of a
        float for either value then that bound is autoscaled.
        """
        inp = inp.split() if type(inp) is str else inp
        self._prop_xrange = "{} {}".format(*inp)
        if any([i == '*' for i in inp]):
            self.ax.autoscale(enable=True, axis='x')
        if inp[0] != '*':
            self.ax.set_xlim(xmin=float(inp[0]))
        if inp[1] != '*':
            self.ax.set_xlim(xmax=float(inp[1]))
    def _set_yrange(self, inp):
        """Set the y range of axis. If an asteriks is detected in place of a
        float for either value then that bound is autoscaled.
        """
        inp = inp.split() if type(inp) is str else inp
        self._prop_yrange = "{} {}".format(*inp)
        if any([i == '*' for i in inp]):
            self.ax.autoscale(enable=True, axis='y')
        if inp[0] != '*':
            self.ax.set_ylim(ymin=float(inp[0]))
        if inp[1] != '*':
            self.ax.set_ylim(ymax=float(inp[1]))

class GraphView(_View):
    """Simple x-y data View class.
    
    The GraphView class is a simple View class for plotting x-y data points. It
    consists of a single axes, and can contain any number of xyData Plot objects.
    It is also the default View mode opened when the program starts.
    """
    def __init__(self, fig):
        """Add axes to the provided figure and set the view mode to 'graph'. The
        available arguments for the user in view mode are defined below, and each
        has a '_set_<argument>' method further down that is called when the
        argument is parsed.
        """
        _View.__init__(self, figure=fig, mode='graph')
        self.ax = fig.add_subplot(111)
        # Line arguments
        self._add_arg('-lw', '--linewidth',  nargs='+', type=float, metavar='float', 
            help="Width of plot line.")
        self._add_arg('-lc', '--linecolour', nargs='+', type=str,   metavar='str',
            help="line colours for plots in figure (any valid matplotlib colour, e.g 'red', 'blue')")
        self._add_arg('-ls', '--linestyle',  nargs='+', type=str,   metavar='str',
            help="line styles for plots in figure (any valid matplotlib style, e.g '-', '--', 'none')")
        self._add_arg('-m',  '--marker',     nargs='+', type=str,   metavar='str',
            help="marker styles for plots in figure (any valid matplotlib style, e.g 'o', 'x', 'none')")
        self._add_arg('-ms', '--markersize', nargs='+', type=str,   metavar='float',
            help="marker sizes for plots in figure (any valid matplotlib style, e.g '-', '--', 'none')")
        self._add_arg('-l',  '--label',      nargs='+', type=str,   metavar='str',
            help="labels for plots in figure legend, use '#' for spaces")
        # Note: argparse does not always handle more exotic variable types and
        # uses in the expected manner. As a result we prefer to use either float
        # or int values, and parse arrays, tuples etc. as multiple arguments, or
        # booleans as strings.

    # View Plot objects
    class Plot(_View.Plot):
        """The Plot object for the GraphView class. Called by cmdPrompt to
        initialise a plot object in the current View instance for a new data file.
        
        """
        def __init__(self, data, ax):
            """Add plot points from Data object to provided axes from the parent 
            View instance. Call the superclass __init__ to store the associated
            Data object and axes for this Plot.
            
            Arguments specified in the parent View instance's __init__ method
            that act on plots have their '_set_' methods defined here. Every
            '_set_' method should first store the input value in the relevant
            '_prop_' attribute to allow the configuration to be written to file.

            """
            _View.Plot.__init__(self, data, ax)
            x, y = self._Data.dat
            self._plot, = ax.plot(x, y)
        # Line methods
        def _set_linewidth(self, inp):
            """Set line width of plot."""
            self._prop_linewidth = inp
            inp = float(inp)
            plt.setp(self._plot, linewidth=inp)
        def _set_linecolour(self, inp):
            """Set line colour of plot."""
            try:
                inp = tuple(float(i) for i in inp.split(','))
            except:
                pass
            self._prop_linecolour = inp
            plt.setp(self._plot, color=inp)
        def _set_linestyle(self, inp):
            """Set line style of plot."""
            self._prop_linestyle = inp
            plt.setp(self._plot, ls=inp)
        def _set_marker(self, inp):
            """Set marker style of plot."""
            if inp == 'none':
                inp = 'None'
            self._prop_marker = inp
            plt.setp(self._plot, marker=inp)
        def _set_markersize(self, inp):
            """Set size of plot markers."""
            self._prop_markersize = inp
            inp = float(inp)
            plt.setp(self._plot, ms=inp)
        def _set_label(self, inp):
            """Set legend labels of plot."""
            self._prop_label = inp
            if inp == 'none': #line has no label
                plt.setp(self._plot, label='__nolegend__')
            else:
                inp = inp.replace('#', ' ')
                plt.setp(self._plot, label=inp)
            plt.legend()

    def add_plot(self, data):
        _View.add_plot(self, data)
        self._set_xrange('* *')
        self._set_yrange('* *')

class StickView(_View):
    """Simple x-y data View class.
    
    The GraphView class is a simple View class for plotting x-y data points. It
    consists of a single axes, and can contain any number of xyData Plot objects.
    It is also the default View mode opened when the program starts.
    """
    def __init__(self, fig):
        """Add axes to the provided figure and set the view mode to 'graph'. The
        available arguments for the user in view mode are defined below, and each
        has a '_set_<argument>' method further down that is called when the
        argument is parsed.
        
        """
        _View.__init__(self, figure=fig, mode='stick')
        self.ax = fig.add_subplot(111)
        #self.ax.set_yscale('log')
        # Line arguments
        self._add_arg('-lw', '--linewidth', nargs='+', type=str, metavar='float',
            help="Width of plot line.")
        self._add_arg('-lc', '--linecolour', nargs='+', type=str,   metavar='str',
            help="line colours for plots in figure (any valid matplotlib colour, e.g 'red', 'blue')")
        self._add_arg('-ls', '--linestyle',  nargs='+', type=str,   metavar='str',
            help="line styles for plots in figure (any valid matplotlib style, e.g '-', '--', 'none')")
        self._add_arg('-m',  '--marker',     nargs='+', type=str,   metavar='str',
            help="marker styles for plots in figure (any valid matplotlib style, e.g 'o', 'x', 'none')")
        self._add_arg('-ms', '--markersize', nargs='+', type=str,   metavar='float',
            help="marker sizes for plots in figure (any valid matplotlib style, e.g '-', '--', 'none')")
        self._add_arg('-l', '--label', nargs='+', type=str,   metavar='str',
            help="labels for plots in figure legend, use '#' for spaces")
        # Note: argparse does not always handle more exotic variable types and
        # uses in the expected manner. As a result we prefer to use either float
        # or int values, and parse arrays, tuples etc. as multiple arguments, or
        # booleans as strings.
    
    # View Plot objects
    class Plot(_View.Plot):
        """The Plot object for the GraphView class. Called by cmdPrompt to
        initialise a plot object in the current View instance for a new data file.
        """
        def __init__(self, data, ax):
            """Add plot points from Data object to provided axes from the parent 
            View instance. Call the superclass __init__ to store the associated
            Data object and axes for this Plot.
            
            Arguments specified in the parent View instance's __init__ method
            that act on plots have their '_set_' methods defined here. Every
            '_set_' method should first store the input value in the relevant
            '_prop_' attribute to allow the configuration to be written to file.
            """
            _View.Plot.__init__(self, data, ax)
            sticks = self._Data.dat
            self._plot = ax.add_collection(
                matplotlib.collections.LineCollection(sticks))
            self._prop_linecolour = self._plot._edgecolors[0]
            self._markers = None
        # Line methods
        def _set_linewidth(self, inp):
            """Set line width of plot."""
            self._prop_linewidth = inp
            if inp == 'none':
                inp = '0'
            inp = float(inp)
            plt.setp(self._plot, linewidth=inp)
        def _set_linecolour(self, inp):
            """Set line colour of plot."""
            try:
                inp = tuple(float(i) for i in inp.split(','))
            except:
                pass
            self._prop_linecolour = inp
            plt.setp(self._plot, color=inp)
        def _set_linestyle(self, inp):
            """Set line style of plot."""
            self._prop_linestyle = inp
            if inp == 'none':
                self._set_linewidth('none')
                return
            elif hasattr(self, '_prop_linewidth'):
                if self._prop_linewidth == 'none':
                    self._set_linewidth(1)
            plt.setp(self._plot, ls=inp)
        def _set_marker(self, inp):
            """Set marker style of plot."""
            self._prop_marker = inp
            if inp == 'none':
                if self._markers:
                    self._markers.remove()
                    self._markers = None
                    return
                else:
                    return
            elif not self._markers:
                self._add_stick_markers()
            plt.setp(self._markers, marker=inp)
        def _set_markersize(self, inp):
            """Set size of plot markers."""
            self._prop_markersize = inp
            inp = float(inp)
            plt.setp(self._markers, ms=inp)
        def _set_label(self, inp):
            """Set legend labels of plot."""
            self._prop_label = inp
            if inp == 'none': #line has no label
                plt.setp(self._plot, label='__nolegend__')
            else:
                inp = inp.replace('#', ' ')
                if self._markers is None:
                    plt.setp(self._plot, label=inp)
                else:
                    plt.setp(self._markers, label=inp)
            plt.legend()

        def _add_stick_markers(self):
            "Plots scatter graph as markers for stick view"
            self._markers, = plt.plot(self._Data.dat[:,1,0], self._Data.dat[:,1,1], 
                                     linestyle='none', linewidth=0,
                                     color=self._prop_linecolour)

class linelistComparisonView(_View):
    def __init__(self, fig):
        _View.__init__(self, figure=fig, mode='resids')
        self.ax = fig.add_subplot(111)
        #Line arguments
        self._add_arg('-lc', '--linecolour', nargs='+', type=str,   metavar='str',
            help="line colours for plots in figure (any valid matplotlib colour, e.g 'red', 'blue')")
        self._add_arg('-m',  '--marker',     nargs='+', type=str,   metavar='str',
            help="marker styles for plots in figure (any valid matplotlib style, e.g 'o', 'x', 'none')")
        self._add_arg('-ms', '--markersize', nargs='+', type=str,   metavar='float',
            help="marker sizes for plots in figure (any valid matplotlib style, e.g '-', '--', 'none')")
        self._add_arg('-l', '--label', nargs='+', type=str,   metavar='str',
            help="labels for plots in figure legend, use '#' for spaces")
    
    class Plot(_View.Plot):
        def __init__(self, data, ax):
            _View.Plot.__init__(self, data, ax) #store data and ax locally as self.data & self.ax
            
            # Ask for second linelist file to compare to
            secFile       = input('Linelist file to compare to: ')
            fType         = detect_filetype(secFile)
            self._secData = fType(secFile).read_file()
            
            # Convert data to pandas dataframe for easy merging
            mergers = [
                'rotational_final', 'rotational_initial', 
                'vibrational_final', 'vibrational_initial'
                ] #column headers to merge on
            self._dfPrim = pd.DataFrame(
                data=self._Data.dat, columns=self._Data.cols
                ).astype(self._Data.typedict)
            self._dfSec = pd.DataFrame(
                data=self._secData.dat, columns=self._secData.cols
                ).astype(self._secData.typedict)
            self._dfPrim.insert(0, 'line_number', 
                np.arange(1, self._dfPrim.shape[0]+1)
                ) #add column with index of line in original file
            self._dfSec.insert(0, 'line_number', 
                np.arange(1, self._dfSec.shape[0]+1)
                )
            # Merge on shortest linelist to minimise surplus lines
            short = 'right' if self._dfSec.shape[0] < self._dfPrim.shape[0] else 'left'
            self.dfComp = self._dfPrim.merge(self._dfSec, on=mergers, how=short)
            # Print merge details
            print("{0} has {1} entries; {2} has {3}. Matched {4} entries"
                .format(
                self._Data.fName, self._dfPrim.shape[0], secFile, 
                self._dfSec.shape[0], self.dfComp.shape[0]
                ))
            self._set_compareon(
                input("Merged, select quantity to compare on: ")
            )
        
        def _set_compareon(self, inp):
            """Select quantity to compare lines on."""
            if hasattr(self, '_plot'):
                self._plot.remove()
            if inp in ['A', 'einstein']:
                #pd.set_option('display.max_rows', 50000, 'display.max_columns', 100, 'display.width', 1000)
                #fOut = open('garbage.txt', 'w')
                #print(self.dfComp, file=fOut)
                yVals = self.dfComp['einstein_A_x']/self.dfComp['einstein_A_y']
            elif inp in ['nu', 'wavenumber']:
                yVals = self.dfComp['wavenumber_x'] - self.dfComp['wavenumber_y']
            else:
                print("Comparison not recognised.")
            xVals = self.dfComp['energy_final_cm_x']
            self._plot, = plt.plot(xVals, yVals, ls='none')

        def _set_linecolour(self, inp):
            """Set line colour of plot."""
            try:
                print(inp)
                inp = tuple(float(i) for i in inp.split(','))
                print(inp)
            except:
                pass
            self._prop_linecolour = inp
            plt.setp(self._plot, color=inp)
        def _set_marker(self, inp):
            """Set marker style of plot."""
            if inp == 'none':
                inp = 'None'
            self._prop_marker = inp
            plt.setp(self._plot, marker=inp)
        def _set_markersize(self, inp):
            """Set size of plot markers."""
            self._prop_markersize = inp
            inp = float(inp)
            plt.setp(self._plot, ms=inp)


