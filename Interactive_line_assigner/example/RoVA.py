################################################################################
# Author: Ryan Brady | Version 1.0.0 | Last Update: 10/09/2022                 #
#                                                                              #
# ABOUT: This code is a line assignment code that allows one to perform the    #
#        full assignment project in one interactive window. This includes      #
#        spectral line profile fitting, quantum number assignment, and an      #
#        inbuilt Combination Difference (CD) test.                             #
# ---------------------------------------------------------------------------- #
# INPUTS: experimental spectrum, theoretical line list, spectrum type params   #
#         in order of command line arguments:                                  #
#                                                                              #
# (1) Experimental spectrum: must have first 2 columns being position/         #
#     intensity.                                                               #
# (2) Exp spectra type: 'high' when line profiles are resolved (high res);     #
#     'stick' when only position and intensity of lines are known.             #
# (3) theoretical line list: must have first 2 columns being position/         #
#     intensity.                                                               #
# (4) theoretical spectrum type: 'high' or 'stick' (see above). When 'high' is #
#     given, both a broadened theoretical spectrum AND a line list should be   #
#     given (within argument 5).                                               #
# (5) Line list to be given when providing aswell a broadened theoretical      #
#     spectrum.                                                                #
#                                                                              #
#  Alternatively, the inputs can be provided within a configuration file       #
#  called 'config.txt'.                                                        # 
#                                                                              # 
# OUTPUTS: 'assignments.txt' which gives line (profile) parameters and QN      #
#          assignments; 'CD.txt' containing all CD test results.               #                                                  
#                                                                              #
################################################################################
#
## ---------------------------- LIBRARY IMPORT: --------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import Button
from lmfit.models import GaussianModel, VoigtModel, LinearModel, ConstantModel
from scipy.special import wofz
import sys
import os
from os.path import exists
# %matplotlib qt (for windows users)
#
## ------------------------ Command line arguments: ---------------------------

if len(sys.argv)==1:
    if os.path.exists('config.txt')==False:
        print(f"'config.txt' doesn't exist...")
        print("Please: (1) provide 'config.txt' in the working directory;")
        print(" (2) provide command line arguments.")
        sys.exit('Now exiting...')

    with open('config.txt','r') as f:
        exp = [line.split()[1] for line in f if "exp-datafile" in line][0]
    with open('config.txt','r') as f:
        exp_res = [line.split()[1] for line in f if "exp-type" in line][0]
    with open('config.txt','r') as f:
        theory = [line.split()[1] for line in f if "calc-datafile" in line][0]
    with open('config.txt','r') as f:
        theory_res = [line.split()[1] for line in f if "calc-type" in line][0]
    with open('config.txt','r') as f:
        LL = [line for line in f if "Line-List" in line]   
    with open('config.txt','r') as f:
        CDf = [line for line in f if "CD" in line]

    if LL != []:
        Line_list = LL[0].split()[1]
    if CDf != []:
        CD = CDf[0].split()[1]

else:
    exp = str(sys.argv[1])
    exp_res = str(sys.argv[2])
    theory = str(sys.argv[3])
    theory_res = str(sys.argv[4])

    if len(sys.argv)==6:
        Line_List = str(sys.argv[5])
    if len(sys.argv)==7:
        CD = str(sys.argv[6])

# Error messages
if os.path.exists(exp)==False:
    print(f"Your input experimental data file '{exp}' can't be found!")
    print()
    print("Please provide a valid file directory. Now exiting...")
    sys.exit()

if os.path.exists(theory)==False:
    print(f"Your input line list data file '{theory}' can't be found!")
    print()
    print("Please provide a valid file directory. Now exiting...")
    sys.exit()   

if (exp_res!="stick") and (exp_res!="high"):
    print("Oops, you didn't specify an accepted experimental data type!")
    print()
    print("Please use the following data types:")
    print("(1) 'stick': only positions and intensity of single lines are known")
    print("(2)  'high': Spectral line profiles are resolvable (broadened)")
    print()
    print("Now exiting...")
    sys.exit()

if (theory_res!="stick") and (theory_res!="high"):
    print("Oops, you didn't specify an accepted theory data type!")
    print()
    print("Please use the following data types:")
    print("(1) 'stick': only positions and intensity of single lines are known")
    print("(2)  'high': Spectral line profiles are resolvable (broadened)")
    print()
    print("Now exiting...")
    sys.exit()

## ------------------------- FUNCTION DEFINITIONS ------------------------------
def load_data(exp_fname, model_fname):                 ## Data loading function
    ''' Data Loading function with 2 arguments.
     ---------------------------------------------------------------------------
    - exp_fname = experimental spectra file name: must be in a space delimited  
                                                  format with energy/wavenumber/
                                                  frequency on first column and
                                                  intensity in second column.

    - model_fname line list file name:  must be in a space delimited format with 
                                        energy/wavenumber/frequency on first 
                                        column and intensity in second column 
                                        followed by QNs.                         
     ---------------------------------------------------------------------------
     The output is a series of global variables:
     - experimental spectra position (x) and intensity (y),
     - line list/model spectra position (nu) and intensity (I),
     - The reduced line list dataframe, concatonated to the region of the exp 
       dataset (df_predicted),
     - The experimental spectra dataframe (df_exp) and column header 
       (exp_header).
     '''
    ## Global variables to be set
    global x                                       # experimental line position
    global y                                       # experimental line intensity
    global nu                                      # theoretical line position
    global I                                       # theoretical line intensity
    global df_predicted                            # line list dataframe
    global df_exp                                  # experimental dataframe
    global exp_header                              # column headers for ^^

    ## Load experimental data                   
    df_exp = pd.read_csv(exp_fname)  # experimental dataframe
    exp_header = df_exp.columns                       # exp data column headers
    x = np.array(df_exp[exp_header[0]])               # experimental energies
    y = np.array(df_exp[exp_header[1]])               # experimental intensities
    
    ## Min and max energies from experiment
    xmin = np.min(x) 
    xmax = np.max(x)
    
    ## Load theoretical data         
    df_mod = pd.read_csv(model_fname)                      # line list dataframe
    mod_header = df_mod.columns                  # line list data column headers
    
    ## truncate the model dataframe by xmin and xmax using boolean masks
    df_predicted = df_mod[ (df_mod[mod_header[0]] <= xmax) 
                         & (df_mod[mod_header[0]] >= xmin) 
                         ]

    #set model spectra postion nu, and intensities I of the truncated dataframe
    nu = np.array(df_predicted[mod_header[0]])
    I  = np.array(df_predicted[mod_header[1]])


def range_finder(n, minn, maxn):                       ## Range finding function
    '''Function with 3 arguments.
     ---------------------------------------------------------------------------
    - n: array of length len(n),
    - minn: minimum value of desired range to be found in array n,   
    - maxn: maximum value of desired range to be found in array n.   
     ---------------------------------------------------------------------------
     The output are the index's of the elements in n closest/containing range 
     values minn and maxn:
     - index[0]: index of element in array n 
     - line list/model spectra position (nu) and intensity (I).
     '''
    # return the indexes where boolean_array = TRUE (i.e. within the range)
    boolean_array = np.logical_and(n >= minn, n <= maxn)
    index = np.where(boolean_array)[0]
    return index[0], index[-1]

                            
def Voigt_amp(p,y):      ## Amplitude of Voigt based on height of spectral line
    '''Function with 2 arguments.
     ---------------------------------------------------------------------------
    - p: parameter object of voigt model to fit to spectral lines,
    - y: array of intensities of spectral line     
     ---------------------------------------------------------------------------
     The output is the corresponding amplitude of the Voigt function given a 
     height np.max(y).
     '''
    z = (complex(0,p['gamma'].value))/(p['sigma'].value*np.sqrt(2))
    return np.max(y)*p['sigma'].value*np.sqrt(2*np.pi)/wofz(z).real
#   
# ----------------------- The 'onselect' functions -----------------------------
#
def onselect_QN(xmin, xmax):              ## Extract QN of line within selection
    '''Onselect function (2 arguments).
     ---------------------------------------------------------------------------
    - xmin: minimum x position of slsection box on interactive plot.
    - xmax: maximum x position of selection box on interactive plot.     
     ---------------------------------------------------------------------------
     The output is a global variable called QN which contains the quantum
     numbers of the line that was 
     selected.
     '''
    # Global variable to be set
    global QN  # Quantum number extracted by line within selection in line list
    
    # Find indices of model line positions within selection xmin and xmax using
    start, end = range_finder(nu, xmin,xmax)

    # extract quantum numbers of row start from model dataframe 
    QN = df_predicted.iloc[[start]]
    
def onselect_line(xmin, xmax, function):   # Extract pos and inten of exp STICK
    '''Onselect function (2 arguments).
     ---------------------------------------------------------------------------
    - xmin: minimum x position of slsection box on interactive plot.
    - xmax: maximum x position of selection box on interactive plot.     
     ---------------------------------------------------------------------------
     The output is a global variable called exp_line_params which contains the 
     experimental position and intensity the experimental line that was selected
     '''
    # Global variable to be set
    global exp_line_params

    # Find indices of experimental line positions within selection xmin and xmax
    start, end = range_finder(x, xmin,xmax)

    # extract exp line position and intensity of row start from exp dataframe
    exp_line_params = df_exp.iloc[[start]]
    
def onselect_line_profile(xmin, xmax):   #Extract voigt line profile parameters 
    '''Onselect function (2 arguments).
     ---------------------------------------------------------------------------
    - xmin: minimum x position of slsection box on interactive plot.
    - xmax: maximum x position of selection box on interactive plot.     
     ---------------------------------------------------------------------------
     This function fits a {voigt + constant} model (equation, see 
     https://en.wikipedia.org/wiki/Voigt_profile) to the data within the limits 
     of the selected x position {xmin, xmax} using the Lmfit library. The output 
     are the global variables x_fit (fitted x-range), y_fit (fitted y-range), 
     fit_out (the fit parameters and results), exp_line_params (the voigt line
     profile parameter = position, height, FWHM, ...), and param_header (the 
     names of the parameters).
     '''

    global x_fit, y_fit, fit_out, exp_line_params, param_header, residual, res_amp

    ## define arrays for the x (pos) and y (inten) to fit with a voigt function
    y_fit = np.array( y[(x < xmax) & (x > xmin)] )
    x_fit = np.array( x[(x < xmax) & (x > xmin)] )
    
    ## Define intitial voigt model and guess start parameters
    model0 = VoigtModel()
    p0 = model0.guess(y_fit, x=x_fit)
    
    ## initialise parameters, setting hard limits for the amplitude
    mod = VoigtModel()+ConstantModel()
    pars = mod.make_params()
    pars.add('amplitude', value = Voigt_amp(p0,y_fit),
              min = 0.95*Voigt_amp(p0,y_fit)
            ) # See note A below
    pars.add('center', value = p0['center'])     # line center position
    pars.add('sigma', value = p0['sigma'])       # Gaussian comp. FWHM
    pars.add('gamma', value = p0['gamma'])       # Lorentzian comp. FWHM
    pars.add('fwhm', value = p0['fwhm'])         # FWHM
    pars.add('c', value = 0)                     # constant shift to line base
    
    ### A: This uses a guess for the amplitude based on the height of the line,
    ### given some buffer (5% here in this case, change this to something that
    ### works better (order of the exp uncertainty)).

    # allow gamma (Lorentzian component of FWHM) to vary
    pars['gamma'].set(value = p0['gamma'].value,vary=True, expr='') 
    
    ## Fit the data with model 'mod'
    fit_out = mod.fit(y_fit, pars, x=x_fit)
    
    # parameter names
    param_header = list(fit_out.params.valuesdict().keys())
    
    # extract best fit line profile parameters

    exp_line_params = np.array(fit_out.params)

    # COmpute the residuals of the fit and its max O-C
    residual = np.array(y_fit-fit_out.best_fit)
    res_amp = np.max(abs(residual))


    
def onselect_zoom(xmin, xmax):       ## zoom in on plot based on user selection 
    ''' On-select function that simultaneously zooms on all subplots.'''

    ## Find maximum y-point (intensity) within selection to scale new plot to
    y_max = np.max(y[(x < xmax) & (x > xmin)])
    
    # Re-plot the selected data with axis limits. 
    ax1.set_xlim(xmin, xmax) # set x axis limits to energy range in selection
    ax1.set_ylim(0, y_max)   # set y axis limits to intensity range in selection
    
    ax2.set_xlim(xmin, xmax) # "
    ax2.set_ylim(0, y_max)   # "
       
    ax3.set_xlim(xmin, xmax) # "
    ax3.set_ylim(y_max,0)    # "
    
    # re-draw/update plot
    fig.canvas.draw_idle()                                  
#
## ------------------------------- DATA LOAD -----------------------------------
load_data(exp,theory)                                 
#
## ------------------------------- PLOTTING ------------------------------------
fig = plt.figure(figsize=(8,6))        # initialise plot figure with size 8x6  "
# set up the first 2 subplots: plots within plots. See the matplotlib manual for
# more info at https://matplotlib.org/2.0.2/users/pyplot_tutorial.html

# plot theory stick
ax3 = fig.add_subplot(313)
line3, = ax3.plot(x,np.zeros(len(x)))
line4 = ax3.stem(nu, I, markerfmt=' ', linefmt='magenta',label="model", 
                     use_line_collection = True
                )
ax3.set_yscale('log')
ax3.set_ylim(0,np.max(I))
ax3.set_xlabel(r'Wavenumbers ($cm^{-1}$)')
ax3.set_ylabel(r'$Log_{10}$ Cross-section ($cm^2$ per molecule)')
ax3.invert_yaxis()
ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, 
               shadow=True, ncol=5, fontsize=8.75
          )


if exp_res == 'stick':
    ax1 = fig.add_subplot(311, sharex = ax3)
    line1 = ax1.plot(x, np.zeros(len(x)), color='black',label="Fit Residuals")
    ax1.xaxis.tick_top()
    ax1.set_title('Press left mouse button and drag to zoom')
    ax1.set_ylabel("Residuals, O-C (Intensity Units)")
    plt.legend()
    
    ax2 = fig.add_subplot(312, sharex=ax3)
    line2 = ax2.stem(x, y, markerfmt=' ', linefmt='blue',
                     label="Experimental Spectrum", 
                     use_line_collection = True
                    )
    ax2.set_ylabel("Intensity (arbitrary units)")
    ax2.set_xticks([])

if exp_res == 'high':
    ax1 = fig.add_subplot(311, sharex = ax3)
    line1 = ax1.plot(x, y, color='black',label="Fit Residuals", alpha=0.3)
    ax1.xaxis.tick_top()
    ax1.set_title('Press left mouse button and drag to zoom')
    ax1.set_ylabel("Residuals, O-C (Intensity Units)")
    plt.legend()
    
    ax2 = fig.add_subplot(312, sharex=ax3)
    line2 = ax2.plot(x, y, color='blue',label="Experimental Spectrum")
    ax2.set_ylabel("Intensity (arbitrary units)")
    plt.setp(ax2.get_xticklabels(), visible=False)



# if theory_res=="stick":
#     ax3 = fig.add_subplot(313)
#     line3, = ax3.plot(x,np.zeros(len(x)))
#     line4 = ax3.stem(nu, I, markerfmt=' ', linefmt='magenta',label="model", 
#                      use_line_collection = True
#                     )
#     ax3.set_yscale('log')
#     ax3.invert_yaxis()
#     ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, 
#                shadow=True, ncol=5, fontsize=8.75
#               )

# if theory_res=="high":
#     ax3 = fig.add_subplot(313)
#     line3, = ax3.plot(x,np.zeros(len(x)))
#     line4 = ax3.plot(nu, I, color='magenta',label="model")
#     ax3.set_yscale('log')
#     ax3.invert_yaxis()
#     ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), fancybox=True, 
#                shadow=True, ncol=5, fontsize=8.75
#               )

# make plots have no whitespace between them and so share axes
plt.subplots_adjust(hspace=0.0)

## Below makes use of the on-select functions... updates ax, ax2, and ax3. 
span1 = SpanSelector(ax2, onselect_line_profile, 'horizontal', useblit=True, 
                     rectprops=dict(alpha=0.5, facecolor='green')
                     )
span3 = SpanSelector(ax3, onselect_QN, 'horizontal', useblit=True, 
                     rectprops=dict(alpha=0.5, facecolor='blue')
                     )
## Buttons!
# button positions
fitax = plt.axes([0.7, 0.025, 0.1, 0.04])    
assignax = plt.axes([0.8, 0.025, 0.1, 0.04]) 

# define buttons and their names, colour e.t.c.
fit_button = Button(fitax, 'Fit', hovercolor='0.975')
assign_button = Button(assignax, 'Assign', hovercolor='0.975')

## these functions will trigger on each event of clicking the button!
# plot the fitted voigt line profile superimposed on data
def fit(event,line_fit = ax2.plot(x,np.zeros(len(x)),'r-'), 
        line_res = ax1.plot(x,np.zeros(len(x)),'r-')):

    if exp_res == 'high':  
        ax2.set_ylim(0,1.2*np.max(y_fit))
        line_fit[0].set_xdata(x_fit)
        line_fit[0].set_ydata(fit_out.best_fit)

        line_res[0].set_xdata(x_fit)
        line_res[0].set_ydata(residual)
        ax1.set_ylim(-1.2*res_amp,1.2*res_amp)

        chi2 = np.sum(residual**2)
        #line_res.set_label(f'Current $\chi^2$={chi2}')
        #ax1.legend()
        fig.canvas.draw_idle()         
    if exp_res == 'stick':
        print("You are trying to fit a profile to a stick...")
        print("please change input type to high-res spectrum.")
fit_button.on_clicked(fit)

# extract line information (profile parameters, quantum numbers) and write them 
# to an assignment file
def assignment(event):
    assign = []                             # initialise empty assignment array
    if exp_res == 'stick':
        for i in np.array(exp_line_params)[0]:
            assign.append(i)    

    if exp_res == 'high':
        ax1.autoscale(enable=True,axis='y')
        ax2.autoscale(enable=False)
        ax3.autoscale(enable=False)
        ax1.plot(x_fit,residual, 'g-')
        ax2.plot(x_fit,fit_out.best_fit, 'g-')
        for i in np.array(exp_line_params):
            assign.append(i)                   
        
    assign.append('-> Assignment ->')

    for i in np.array(QN)[0]:
        assign.append(i)                     
         
    ## write the assignment to the file 'assignments.txt'
    with open('assignments.txt', 'a') as f:
        f.write(f"{assign}")
        f.write("\n")
assign_button.on_clicked(assignment)

plt.show()