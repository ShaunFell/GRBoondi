## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

import matplotlib.pyplot as plt
import os, configparser, sys
import numpy as np
import pandas as pd

from Source.OneD import *
from Source.Common.Utils import get_config

## Load config file
config = get_config(sys.argv) #includes error checking

#if the plot and movie paths dont exist, create them
create_output_dirs(config)

#Instantiate printer function that accepts a verbosity level
verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))
verbPrint("Verbosity: ", verbPrint.verbosity)


def main():
    """
    Main function to plot the data

    Raises:
        IOError: Data path not found
        IOError: Data file not found
        RuntimeError: No plot variables specified
        RuntimeError: Number of linestyles does not match number of plot variables.
        RuntimeError: Number of xlims does not match number of plot variables.
    """

    #Extract basic file information
    datapath = config["Header"].get("data_path", fallback = "./data")
    filename = config["Header"].get("data_filename", fallback = "Integrals.dat")

    # get the filename that contains the data
    data_filename_abs = os.path.join(datapath, filename)
    verbPrint("Data file: {0}".format(data_filename_abs))

    #check if data directory exists
    if not os.path.exists(datapath):
        raise IOError("Data path not found: {0}".format(datapath))

    # check if file exists, if not raise SystemError
    if not os.path.exists(data_filename_abs):
        raise IOError("Data file not found: {0}".format(data_filename_abs))
    
    # get the file object
    data_obj = get_data_dataframe(data_filename_abs)

    #get the plot variables
    plot_variables = config["VariableData"].get("plot_variables", "").split()
    if not plot_variables:
        raise RuntimeError("No plot variables specified.")  
 
    
    verbPrint("Plotting variables: {0}".format(plot_variables))

    #get the line styles
    linestyles = config["VariableData"].get("linestyles", "").split()
    if linestyles == "":
        linestyles = None

    if len(linestyles) != len(plot_variables) and linestyles is not None:
        raise RuntimeError("Number of linestyles does not match number of plot variables.")

    # get the x limits
    xlims = np.float64(config["VariableData"].get("xlim", "").split())
    set_xlims = config["VariableData"].getboolean("set_xlim", False)
    if not len(xlims) == len(plot_variables) and set_xlims:
        raise RuntimeError("Number of xlims does not match number of plot variables.")
    if not set_xlims:
        xlims = None

    #setup the style of the pyplot plot
    setup_pyplot(config)
    
    #generate all the plots and save to disk
    make_plots(config, data_obj, plot_variables, xlims, linestyles)
    

if __name__ == "__main__":
    main()
elif __visit_script_file__ == __visit_source_file__:
    raise RuntimeError("This file should not be run with VisIt, only Python!")
