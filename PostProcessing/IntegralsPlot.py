import matplotlib.pyplot as plt
import os, configparser, sys
import numpy as np
import pandas as pd


#Setup the config file parser
config = configparser.ConfigParser()
config.read(sys.argv[1])

#Extract basic file information
plotpath = config["Output"]["output_plot_path"]
datapath = config["Header"]["data_path"]
filename = config["Header"]["data_filename"]
verbosity = config["Header"].getint("verbosity",0)

#If the plot path doesnt exist, create it
if not os.path.exists(plotpath):
    os.mkdir(plotpath)

#Setup print function that respects user-defined verbosity
def verbPrint(*objects):
	if verbosity:
		print(*objects)

verbPrint("Output plots directory: {0}".format(plotpath))


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
    #import the data to a pandas dataframe
    return pd.read_csv(filepath, delim_whitespace=True, names=header, engine="python", header=0, index_col=False)

def setup_plots():
    """
    Setup the plot style
    """
    # turn off interactive plotting and set plot style
    plt.ioff()
    plt.style.use(config["PlotConfig"].get("plotstyle", "default"))

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

def save_fig(path, varname):
    """
    Save the pyplot figure

    Args:
        path (str): absolute path to the output directory
        varname (str): name of the variable
    """
    plt.savefig(os.path.join(path,  varname), format = config["Output"].get("fileform", "png").lower(), dpi = config["Output"].getint("dpi", 600))

def setup_figure(): 
    """
    Setup the figure environment
    """
    plt.figure()
    plt.tight_layout()


def make_plot(data_obj, plot_variable, plotbounds = None, linestyle = '-' ):
    """
    Make the actual plot

    Args:
        data_obj (pd.DataFrame): a pandas dataframe holding the data
        plot_variable (str): the name of the variable to be plotted
        plotbounds (list, optional): the plot bounds. Defaults to None.
        linestyle (str, optional): the plot line style. Defaults to '-'.
    """

    #get the scaling
    plotscaling = config["PlotConfig"].get("xscaling", "linear") + config["PlotConfig"].get("yscaling", "linear")

    #get the plot function as determined by the scaling and plot it
    select_plot_func(plotscaling)(data_obj["time"], data_obj[plot_variable], linestyle = linestyle, label=plot_variable)

    #set the labels
    plt.xlabel(config["PlotConfig"].get("xlabel", "Time"))
    plt.ylabel(plot_variable)

    #if the plot bounds are specified, set them
    if plotbounds is not None:
        plt.xlim(plotbounds[0], plotbounds[1])


def make_plots(data_obj, plot_variables, plotbounds = None,  linestyles = None):
    """
    Make multiple plots

    Args:
        data_obj (pd.DataFrame): a pandas dataframe holding the data
        plot_variables (list): the names of the variables to be plotted, in a list
        plotbounds (list, optional): the plot bounds. Defaults to None.
        linestyles (list, optional): the plot line styles in a list format. Defaults to None.
    """

    #flag to tell if the plots should all be in one figure
    oneplot = config["VariableData"].getboolean("one_plot", False)
    if not oneplot:

        #iterate over the variables and plot them using the single plot function
        for i in range(len(plot_variables)):
            setup_figure()
            plotvar = plot_variables[i]

            if linestyles is not None:
                linestyle = linestyles[i]
            else:
                linestyle = '-'

            make_plot(data_obj, plotvar, plotbounds, linestyle)
            plt.legend()
            save_fig(plotpath, plotvar)

    #combine all the plots into one
    else:
        setup_figure()
        for i in range(len(plot_variables)):
            plotvar = plot_variables[i]

            if linestyles is not None:
                linestyle = linestyles[i]
            else:
                linestyle = '-'
            make_plot(data_obj, plotvar, plotbounds, linestyle)
        plt.legend()
        plt.ylabel("") #remove ylabel for combined plot
        save_fig(plotpath, "combined")


def main():
    """
    Main function to plot the data

    Raises:
        SystemError: Data file not found
        SystemError: No plot variables specified
        SystemError: Number of x limits does not match number of plot variables
        SystemError: Number of linestyles does not match number of plot variables
    """

    # get the filename that contains the data
    data_filename_abs = os.path.join(datapath, filename)
    verbPrint("Data file: {0}".format(data_filename_abs))

    # check if file exists, if not raise SystemError
    if not os.path.exists(data_filename_abs):
        raise SystemError("Data file not found: {0}".format(data_filename_abs))
    
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
        raise SystemError("Number of linestyles does not match number of plot variables.")

    # get the x limits
    xlims = np.float64(config["VariableData"].get("xlim", "").split())
    set_xlims = config["VariableData"].getboolean("set_xlim", False)
    if not len(xlims) == len(plot_variables) and set_xlims:
        raise SystemError("Number of xlims does not match number of plot variables.")
    if not set_xlims:
        xlims = None

    #setup the plot style
    setup_plots()
    
    #make the plots
    make_plots(data_obj, plot_variables, xlims, linestyles)

    




if __name__ == "__main__":
    main()