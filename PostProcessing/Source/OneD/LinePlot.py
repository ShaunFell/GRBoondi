## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

from Source.Common.Utils import *
try:
    import matplotlib.pyplot as plt
except:
    raise RuntimeError("matplotlib not found. Note: this file should not be loaded by VisIt.")

def make_plot(config, data_obj, plot_variable, plotbounds = None, linestyle = '-' ):
    """
    Make a line plot, given pandas dataframe and plot variable

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

def make_plots(config, data_obj, plot_variables, plotbounds = None,  linestyles = None):
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
            setup_pyplot_figure()
            plotvar = plot_variables[i]

            if linestyles is not None:
                linestyle = linestyles[i]
            else:
                linestyle = '-'

            #save this single plot to disk
            make_plot(config, data_obj, plotvar, plotbounds, linestyle)
            plt.legend()
            save_pyplot_fig(config, plotvar)

    #combine all the plots into one
    else:
        setup_pyplot_figure()
        
        #iterate over the variables and plot them using the single plot function
        for i in range(len(plot_variables)):
            plotvar = plot_variables[i]

            if linestyles is not None:
                linestyle = linestyles[i]
            else:
                linestyle = '-'
            
            make_plot(config, data_obj, plotvar, plotbounds, linestyle)
        
        #save all the plots in a single figure environment to disk
        plt.legend()
        plt.ylabel("") #remove ylabel for combined plot
        save_pyplot_fig(config, "combined")
        