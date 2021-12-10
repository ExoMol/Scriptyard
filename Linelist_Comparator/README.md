# Linelist Comparator
_by Wilfrid Somogyi_

This directory contains a small python module that can be used for comparing spectral linelists from two different sources according to one or more quantities on a line-by-line basis. Natively supported formats are the Exomol `.states` + `.trans` format, and the Hitran 2004 120 character `.par` format. Any other arbitrary linelist format can be specified by populating the first line of the relevant linelist file with column headers for the relevant quantities and quantum numbers. Any of the following are currently recognised quantities.

## Specifying a linelist

### State Quantities
These are quantities associated with distinguishable quantum states. As such, each line in the linelist file will have two of each quantity associated with it - namely one for the initial state, and one for the final state. In each case the name of the quantitiy given below should be suffixed by either `_i` or `_f` accordingly, i.e label the linelist column corresponding to the initial state energy `energy_i`, the final state energy `energy_f`, and so on.

* `energy`
* `lifetime`
* `degeneracy`
* `angmom_total`
* `angmom_electronic`
* `angmom_orbital`
* `angmom_spin`
* `angmom_proj_total`
* `angmom_proj_electronic`
* `angmom_proj_orbital`
* `angmom_proj_spin`
* `vibrational`
* `parity_total`
* `parity_rotationless`
* `electronic_state`
* `state_number`

### Transition Quantities
These are quantities associated with transition between two quantum states. As such each line in the linelist file will have only one of each quantity, and no additional suffix is required.
* `transition_wavenumber`
* `einstein_coefficient`
* `transition_linestrength`
* `transition_intensity`

## Using `llcomp`

### Loading a linelist
To load a linelist from file(s), use one of the following functions. This will create an `llcomp` object, which allows the user to perform various operations on the data.

* `llcomp.file_to_linelist(fname)`
  - Expects the name of a single linelist file, where each row corresponds to a transition. Requires user-defined columns headers, from the above list of recognised quantities, as the first line of the file.
* `llcomp.exomol_to_linelist(states_file=None, trans_file=None)`
  - Expects a linelist in the two file Exomol format. Does not require user-defined column headers.
* `llcomp.hitran_to_linelist(fname)`
  - Expects a linelist in the Hitran 2004 format. Does not require user-defined column headers.

### Filtering data
To filter data in a `Linelist` object, apply the `filter_data()` method. Multiple filters can be applied simultaneously by providing a list, for example:

```
mylinelist = llcomp.file_to_linelist("myfile.txt")
mylinelist.filter_data([["angmom_electronic_i", "==", 2], ["vibrational", ">", 1]])
```

Note also that filters can be applied to either the initial or final state using the relevant prefix, or to both by writing the label with no prefix (e.g `vibrational` in the example above.

### Comparing linelists
To compare two linelists, one must create a `llcomp.mergedLinelist` instance. This is done by providing the two `Linelist` objects you would like to compare, e.g

```
mylinelist = llcomp.file_to_linelist("myfile.txt")
exomollinelist = llcomp.exomol_to_linelist(states_file="linelist.states", trans_file="linelist.trans")
comparelist = llcomp.MergedLinelist(mylinelist, exomollinelist)
```
By default `llcomp` will merge transitions according to the values of `angmom_total_i`, `angmom_total_f`, `vibrational_i`, `vibrational_f`, `electronic_state_i` and `electronic_state_f`. Remaining quantities will then be appended with `_L` or `_R` depending on whether they belong to the left linelist or the right linelist (`mylinelist` and `exomollinelist`, respectively, in the example above). 

