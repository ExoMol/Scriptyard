# Molpro Extractor
_by Wilfrid Somogyi_

This script is a template for extracting quantities from Molpro output files. The top of the file `extract.py` contains an example setup, with the Molpro state numbers for the energies to be extracted as well as the couplings between them.

The script leverages regex expressions to detect specific sections of the Molpro output, placing the numerical value to be extracted in the first capture group. For example, the regex key

```python
f"""!MRCI STATE\s*{re.escape(state)}\s*Energy\s*(-?\d{{3}}\.\d{{12}})i?"""
``` 

Searches for a line containing `!MRCI STATE` followed by any number of blank spaces, and then the value of the variable `state` followed again by any number of spaces. The regex in brackets `()` is a capture group that matches a positive or negative number with three digits before the decimal place, and 12 afterwards. The script assumes only a single capture group per trigger, corresponding to the value to be extracted.

Currently the script only works with MRCI calculations, but it is trivial to add other calculations simply by including a regex expression to match the relevant line in the Molpro output file.
