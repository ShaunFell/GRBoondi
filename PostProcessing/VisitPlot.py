#
## Run this using 'visit -nowin -cli -s 2DPlot.py paramfile'
## to run in parallel do 'visit -np 32 -nowin -cli -s 2DPlot.py paramfile'
## NOTE: you will need ffmpeg installed to make the movie
#

import glob, os, sys, configparser
import numpy as np
from Source.Common.Utils import VerbosityPrint, MultipleDatabase, Close_Database, get_config
from Source.Common.Utils import create_output_dirs
from Source.Common.Engine import *
from Source.Common.Plotter import make_slice_plots

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
		OSError: ffmpeg failed. Could not make the movie.
	"""

	# get list of all hdf5 plot files
	hdf5files = glob.glob(os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"] + "*.hdf5"))
	hdf5files.sort() # sort the files by number
	print("Running Plot routines!")

	# get list of plot variables
	plot_variables = config["VariableData"].get("plot_variables", "").split()
	if len(plot_variables) == 0:
		raise SystemError("No plot variables specified!")
	
	verbPrint("Number of plot variables: ", len(plot_variables))
	verbPrint("Variables to be plotted: ", config["VariableData"].get("plot_variables", ""))

	# setup the compute engine
	setup_engine(config)

	#iterate over the plot variables, and create plots for each hdf5
	for i in range(len(plot_variables)):
		plotvar = plot_variables[i]
		
		#get bounds for this plot variables
		config_setplotbounds_string = "set_plot_bounds_" + str(i+1)
		config_plotbounds_string = "plot_variables_range_" + str(i+1)
		print("Plotting slices for {0}...".format(plotvar))
		setplotbounds = config["VariableData"].getint(config_setplotbounds_string, 0)
		plotbounds = tuple(np.float64(config["VariableData"].get(config_plotbounds_string, "0 0").split()))

		#create the slice plot
		plottype = config["Header"].get("plot_type", fallback = '2d')
		make_slice_plots(config, plotvar, hdf5files, plot_type = plottype, setplotbounds = setplotbounds, plotbounds = plotbounds)

		#make a movie if asked
		if MultipleDatabase(config) & config["Output"].getint("make_movie", fallback = 0):
			print ("Making a movie...")
			ffmpeg_cmd = 'ffmpeg -r ' + str(config["Output"].get("movie_framerate", fallback = "5")) + ' -v ' + str(config["Output"].get("verbosity", fallback = "0")) + ' -s 1920x1080 -i ' + config["Output"]["output_plot_path"] +"/"+ plotvar + '%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ' + config["Output"]["output_movie_path"] +"/"+ plotvar+ '.mp4'
			system_status = os.system(ffmpeg_cmd)
			
			#Check ffmpeg status to ensure command executed successfully
			if not system_status == 0:
				raise 	OSError("ffmpeg failed. Could not make the movie. cmd: {0}".format(ffmpeg_cmd))

		print("I've finished!")
	
	#cleanup	
	os.remove("./visitlog.py")
	Close_Database(config) #includes error handling

	print("All Done!")
	sys.exit(0)

if __visit_script_file__ == __visit_source_file__:
	main()
