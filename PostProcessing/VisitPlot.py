#
## Run this using 'visit -nowin -cli -s 2DPlot.py paramfile'
## to run in parallel do 'visit -np 32 -nowin -cli -s 2DPlot.py paramfile'
## NOTE: you will need ffmpeg installed to make the movie
#

import glob, os, sys, configparser
import numpy as np
from Source.Common.Utils import VerbosityPrint, MultipleDatabase, Close_Database, get_config
from Source.Common.Utils import create_output_dirs, get_hdf5_file_list, make_movie
from Source.Common.Engine import *
from Source.Common.Plotter import make_visit_plot

from visit import *

## Load config file
config = get_config(sys.argv) #includes error checking

#if the plot and movie paths dont exist, create them
create_output_dirs(config)

#Instantiate printer function that accepts a verbosity level
verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))
verbPrint("Verbosity: ", verbPrint.verbosity)

verbPrint("HDF5 file path: " + config["Header"]["hdf5_path"])
verbPrint("Output Plots directory: " + config["Output"]["output_plot_path"])
verbPrint("Output Movie Directory: " + config["Output"]["output_movie_path"])

def main():
	"""main function to run VisIt plot routines

	Raises:
		RuntimeError: No plot variables specified!
	"""

	# get list of all hdf5 plot files
	hdf5files = get_hdf5_file_list(config)
	print("Running Plot routines!")	

	# get list of plot variables
	plot_variables = config["VariableData"].get("plot_variables", "").split()
	if not plot_variables:
		raise RuntimeError("No plot variables specified!")
	
	verbPrint("Number of plot variables: ", len(plot_variables))
	verbPrint("Variables to be plotted: ", config["VariableData"].get("plot_variables", ""))

	# setup the compute engine
	setup_engine(config)

	#iterate over the plot variables, and create plots for each hdf5
	for i in range(len(plot_variables)):
		plotvar = plot_variables[i]
		
		#get bounds for this plot variable
		config_setplotbounds_string = "set_plot_bounds_" + str(i+1)
		config_plotbounds_string = "plot_variables_range_" + str(i+1)
		setplotbounds = config["VariableData"].getint(config_setplotbounds_string, 0)
		plotbounds = tuple(np.float64(config["VariableData"].get(config_plotbounds_string, "0 0").split()))


		#create the slice plot
		print("Plotting slices for {0}...".format(plotvar))
		plottype = config["Header"].get("plot_type", fallback = '2d')
		make_visit_plot(config, plotvar, hdf5files, plot_type = plottype, setplotbounds = setplotbounds, plotbounds = plotbounds)

		#make a movie if asked
		make_movie(config, plotvar)

		print("I've finished!")
	
	#cleanup	
	os.remove("./visitlog.py")
	Close_Database(config) #includes error handling

	sys.exit(0)

if __visit_script_file__ == __visit_source_file__:
	main()
