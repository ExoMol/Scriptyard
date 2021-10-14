_by Wilfrid Somogyi_

This is a glorified script for plotting stick spectra, I made it to have a fast way of interactively viewing spectra in a stick format, and quickly creating publication quality stick plots.

It uses Python 3, and requires matplotlib, numpy and pandas libraries, and essentially acts as a glorified wrapper for matplotlib that allows you to save plot setups with a simplified configuration file.

To launch the plotter use the Python module syntax:

```
python -m cmdGraph
```

There are two type of commands:
- Console commands: control the program itself, e.g loading `.stick` files, saving configurations and printing graphs.
  * These are simply keywords, the available console commands can be seen by typing 'help' in the prompt.
- Figure commands: control the appearance of the plot, e.g line styles, colours, labels, etc.
  * These are prefixed by a dash (`-`) or double dash (`--`), the available figure commands can be seen by typing '--help'.
 
## Example of Plotting a `.stick` File
To add a stick file to the plot type
```
> adat HITRAN__16O2__QM.stick
```

And set the appropriate x/y ranges
```
\> --xrange 0 2e4
\> --yrange 0 1e-20
```

## Saving and Loading Plots
After customising your plot you can save it to a file using the save command, 
```
\> save example.cmg
```

The save file consists of the series of commands that are required to replicate the plot, when calling the ``load`` command, the program simply executes these commands, as a result you can create the plot in advance by entering the relevant commands in a text file with the header `---cmdGraph---`.

## Warnings
If a command is entered incorrectly, the program is liable to crash as it has no error handling functions. I usually execute the `save` command at regular intervals when constructing a large plot.

Currently the command 'tight' which executes `plt.tight()` cannot be saved. If you load a figure for which this command is required you will have to manually execute `tight` after loading the saved figure.
