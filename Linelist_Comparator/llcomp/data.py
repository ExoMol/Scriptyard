import pandas as pd
import numpy  as np

def y_as_fx(dataframe, x=None, y=None):
    """Return data series x and y from the dataframe as a (n, 2) ndarray,
    sorted by increasing x value."""
    xy = dataframe[[x, y]].to_numpy()
    return xy[xy[:,0].argsort()]

"""
@todo: Add method for merging linelists.
@body: Add a Linelist method for merging self to another Dataframe via this method
"""
def compare_dataframes(left_df, right_df, merge_on):
    """Internal method for retrieving comparisons to other linelists"""
    return left_df.merge(right_df,
        how='inner',
        on=merge_on,
        suffixes=("_L", "_R")
    )

def detect_file_headers(filename, headers_to_detect):
    """Detect the headers in the first line of a file from a given list.
    arguments
        filename : str
            Name of file to sample first line from.
        heads_to_detect : list-like of strings 
            List of header strings to detect.
    returns
        use_these_columns : list of lists
            List of recognised column headers (str) and their indices (int).
        garbage : list of list
            List of unrecognised column headers (str) and their indices (int).
    """
    use_these_columns = [] #collect recognised columns
    garbage = []
    with open(filename, 'r') as f:
        line = f.readline().rstrip('\n')
        for w, word in enumerate(line.split()):
            if word in headers_to_detect: #if recognised word
                use_these_columns.append([word, w]) #keep word and column index
            else:
                garbage.append([word, w])
        if len(use_these_columns) == 0: #if none known
            return None, garbage
        else:
            return use_these_columns, garbage

def is_iterable(obj, strings=False):
    """Check if object is iterable, return boolean result.
    arguments
        obj : object
            Any Python object.
        strings : bool
            If false, strings don't count as iterables.
    """
    try:
        iter(obj)
    except Exception:
        return False
    else:
        if type(obj) is str:
            return False
        else:
            return True

def print_linelist(linelist, fname="blah.txt", cols=None, num_rows=50):
    df = linelist.dataframe
    pd.set_option('display.max_columns', df.shape[1])
    pd.set_option('display.max_rows', df.shape[0])
    pd.set_option('display.width', 1000)
    with open(fname, 'w') as f: 
        if cols:
            print(df[cols], file=f)
        else:
            print(df, file=f)

branch_dict = {
    "O" : -2,
    "P" : -1,
    "Q" : 0,
    "R" : 1,
    "S" : 2
}

def convert_from_branch(quanta_initial, branch):
    """Calculate final state quanta from initial state quanta and branch label.
    arguments
        quanta_initial : float, int
            The quantum number of the initial state.
        branch :  str
            The transition branch ['O, 'P', 'Q', 'R', 'S'].
    returns
        quanta_final : float, int
            The quantum number of the final state.
    """ 
    if branch in branch_dict:
        quanta_final = quanta_initial + branch_dict[branch]
    else:
        quanta_final = None
    return quanta_final
