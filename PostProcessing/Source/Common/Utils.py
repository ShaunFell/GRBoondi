## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

import os, glob
import matplotlib.pylot as plt
import pandas as pd



def PlotFiles(config):
    """	Find all the plot files and return them as list of absolute path strings
    
    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters

    Raises:
        FileNotFoundError: Error raised if no plot files were found

    Returns:
        list: list of plot files
    """

    #path to hdf5 files plus plot file header string
	filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])

    #get list of file using regex
	files = glob.glob(filename_prefix+ "*.3d.hdf5")

    #filter the file list to ensure regex-captured files are actual plot files
	plot_files = [x for x in files if config["Header"]["plot_header"] in x]

    #ensure plot files exist in provided directory
	if len(plot_files) == 0:
		raise FileNotFoundError("No Plot Files Found!")
	
	return plot_files

def MultipleDatabase(config):
    """Determine whether multiple plot files need to be opened

    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters

    Returns:
        boolean: boolean flag determining if multiple plot files need to be opened
    """
	
	if len(PlotFiles(config))>1:
		print("Opening multiple plot files")
		return True
	else:
		print("Opening single plot file")
		return False

def get_data_dataframe(filepath):
    """
    Read in the data file and return a pandas dataframe

    Args:
        filepath (str): absolute path to the data file

    Returns:
        pd.DataFrame: a pandas dataframe
    """

    #get the header lines
    #from GRBoondi, they should always start with a '#'
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("#"):
                header = line
            else:
                break #stop the loop if theres no more header liens
    
    header = header[1:].strip().split()
    verbPrint("Header: {0}".format(header))
    #import the data to a pandas dataframe and return it
    return pd.read_csv(filepath, delim_whitespace=True, names=header, engine="python", header=0, index_col=False)

def setup_pyplot(config):
    """
    Setup the plot style
    """
    # turn off interactive plotting and set plot style
    plt.ioff()
    plt.style.use(config["PlotConfig"].get("plotstyle", "default"))

def setup_pyplot_figure(): 
    """
    Setup the figure environment
    """
    plt.figure()
    plt.tight_layout()

def select_plot_func(scaling):
    """
    Select the plot function based on the scaling

    Args:
        scaling (str): "linear", "loglog", "loglinear", "linearlog"

    Returns:
        function: the plot function from pyplot
    """
    if scaling == "loglog":
        return plt.loglog
    elif scaling == "loglinear":
        return plt.semilogx
    elif scaling == "linearlog":
        return plt.semilogy
    else:
        return plt.plot

def save_pyplot_fig(config, path, varname):
    """
    Save the pyplot figure

    Args:
        path (str): absolute path to the output directory
        varname (str): name of the variable
    """
    plt.savefig(os.path.join(path,  varname), format = config["Output"].get("fileform", "png").lower(), dpi = config["Output"].getint("dpi", 600))



class VerbosityPrint:
    """
    class that allows printing using user-specified verbosity levels
    
    Instantiate class using a given verbosity level, then call the instance 
    by passing objects to be printed as arguments.

    Example:

        verbPrint = VerbosityPrint(1)
        verbPrint("I'm an example!")

        output: "I'm an example!"

        verbPrint = VerbosityPrint(0)
        verbPrint("I'm another example!")

        output: ""
    """

    def __init__(self, verbosity):
        self.verbosity = verbosity

    def __call__(self, *objects):
        if self.verbosity:
            print(*objects)