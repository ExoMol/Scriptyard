# Interactive Line Assigner
_by R. P. Brady_

This is a python script for the assignment of spectral lines in an interactive way through use of the many cool features of the 'matplotlib.widgets' library. The central idea of this script is to assign the brightest lines using a line list determined _a priori_ to the assignment project through observing common structures between the experimental spectrum and the (potentially) un-refined line list. With the brightest lines assigned, one can confirm the assignments using Combination Difference (CD) tests (to be included in the automated pipeline here), to which the confirmed line positions can be used to refine the line list and improve the quality of further assignments of weaker lines.

## Setup

The python script requires the following libraries/modules to be installed [`numpy`,`pandas`,`matplotlib`,`Lmfit`,`scipy`,`os`,`sys`] which can be done through the following pip install via command line:

```
pip install <module name>
```

### Virtual environment (Not necessary for running the line assigner program):
I would reccomend using a virtual environment when using this script to prevent any potential conflicts with other libraries installed (although not neccesary). You can do this through the entering the following terminal command:

#### Mac:

```
conda create –n <some name> python=3
```

#### Windows:

Open an Anaconda Terminal and enter:

```
conda create –n <some name> python=3 libpython m2w64-toolchain
```

This will create a python environment completely clean of any libraries you have installed on your native machine. To install required libraries please see the previous step above.

## The Input: Data set up and keywords

Two files are needed: (1) the experimental datafile including positions and intensities; (2) the theoretical line list including positions, intensities, and quantum numbers. These files need to be in a `.csv` format.

The script also needs input to where the experiment and theory datafiles are as well as some keywords telling the script what data-type you have provided (i.e. stick spectra or resolved spectra). The input parameters are as follows:

* `exp-datafile` = "~/the/directory/to/your/experimental/datafile/exp.csv"
* `exp-type` = high/stick. 'high' means high resolution experimental data (can see resolved line profiles), 'stick' means only positions and intensities of lines are known.
* `calc-datafile` = "~/the/directory/to/your/theoretical/datafile/theory.csv"
* `calc-type` = high/stick. 'high' means you have modelled your line list with a line profile, 'stick' means stick spectra. NOTE: only stick is supported at the moment, I plan on incoorporating profile modelled line lists soon.

We can input the parameters into the script two ways...

### 1: Command line argument

We run the script with the input parameters as follows:
```
python3 RoVa.py exp-datafile exp-type calc-datafile calc-type
```

### 2: Input file

We have an input file called `config.txt` that contains the keywords listed above followed by your desired input. For example:

```
exp-datafile ./data/H2CO_exp_spec_6546-6900_cm-1.csv 
exp-type high 
calc-datafile ./data/H2CO_calc_dataframe_simplified.csv 
calc-type stick 
```

The order of the keywords doesn't matter.

### Technical detail:

For windows users, please go to line 43 and uncomment it. You should set the backend to `%matplotlib qt` since the interactive plotting window won't pop up. You don't have to do this if you are a mac user, but if the code throws a fit then try `%matplotlib widget`.

## Useage:
Once ran, the script should pop up a window that shows you three windows (please see image in the example folder). The top subplot is used for showing residuals of fitted line profile, the middle plot shows the experimental spectrum, and the bottom plot shows the theoretical line list. How to get assignments is done through the following:

* Zoom in on the region of interest using the `matplotlib` zoom tool,
* Drag over the experimental line you wish to assign in the middle plot, a green box should appear and will be the rsange of data used to fit  functional form (such as a `voigt` curve) to,
* Press the `Fit` button, you should now see a red curve overlaying the experimental spectra which is your fitted line profile. In the top plot you will see the `observed - calculated` / `residuals`. You can re-select and fit as many curves as you want untill you are happy with your fit,
* Drag over the theoretical stick line you want to assign to the experimental line (the line list isn't perfect so you will have to observe for common structures to sus out the correct line... start with the brightest lines),
* Press the `Assign` button.

You should now have a file called `assignments.txt` in your working directory. This will contain all of your experimental line profile parameters such as position, height, FWHM, and area followed by the quantum number assignments.

