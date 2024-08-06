## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

import matplotlib.pyplot as plt
import os, configparser, sys
import numpy as np
import pandas as pd

from 1D import *


def main():
    """
    Main function to plot the data

    Raises:
        SystemError: Data file not found
        SystemError: No plot variables specified
        SystemError: Number of x limits does not match number of plot variables
        SystemError: Number of linestyles does not match number of plot variables
    """

    #Setup the config file parser
    config = configparser.ConfigParser()
    config.read(sys.argv[1])

    #Extract basic file information
    plotpath = config["Output"]["output_plot_path"]
    datapath = config["Header"]["data_path"]
    filename = config["Header"]["data_filename"]

    #If the plot path doesnt exist, create it
    if not os.path.exists(plotpath):
        os.mkdir(plotpath)

    #Instantiate printer function that accepts a verbosity level
    verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))
    verbPrint("Verbosity: ", verbPrint.verbosity)

    # get the filename that contains the data
    data_filename_abs = os.path.join(datapath, filename)
    verbPrint("Data file: {0}".format(data_filename_abs))

    # check if file exists, if not raise SystemError
    if not os.path.exists(data_filename_abs):
        raise IOError("Data file not found: {0}".format(data_filename_abs))
    
    # get the file object
    data_obj = get_data_dataframe(data_filename_abs)

    #get the plot variables
    plot_variables = config["VariableData"].get("plot_variables", "").split()    
    
    #check the plot variables are non-empty
    if not plot_variables:
        raise SystemError("No plot variables specified.")
    
    verbPrint("Plotting variables: {0}".format(plot_variables))

    #get the line styles
    linestyles = config["VariableData"].get("linestyles", "").split()
    if linestyles == "":
        linestyles = None

    if not len(linestyles) == len(plot_variables) and linestyles != None:
        raise RuntimeError("Number of linestyles does not match number of plot variables.")

    # get the x limits
    xlims = np.float64(config["VariableData"].get("xlim", "").split())
    set_xlims = config["VariableData"].getboolean("set_xlim", False)
    if not len(xlims) == len(plot_variables) and set_xlims:
        raise RuntimeError("Number of xlims does not match number of plot variables.")
    if not set_xlims:
        xlims = None

    #setup the style of the pyplot plot
    setup_plots()
    
    #generate all the plots and save to disk
    make_plots(data_obj, plot_variables, xlims, linestyles)
    

if __name__ == "__main__":
    main()