#
## Run this using 'visit -nowin -cli -s 2DPlot.py paramfile'
## to run in parallel do 'visit -np 32 -nowin -cli -s 2DPlot.py paramfile'
## NOTE: you will need ffmpeg installed to make the movie
#

import glob, os, sys, configparser, re
import numpy as np

## Load config file
config = configparser.ConfigParser()
config.read(sys.argv[1])


# Extract config variables
plotpath = config["Output"]["output_plot_path"]
moviepath = config["Output"]["output_movie_path"]

#if the plot and movie paths dont exist, create them
if not os.path.exists(plotpath):
	os.mkdir(plotpath)
if not os.path.exists(moviepath):
	os.mkdir(moviepath)

print("HDF5 file path: " + config["Header"]["hdf5_path"])
print("Output Plots directory: " +plotpath)
print("Output Movie Directory: " +moviepath)

#Instantiate printer function that accepts a verbosity level
verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

verbPrint("Verbosity: ", verbPrint.verbosity)

def main():

	# get list of all hdf5 plot files
	hdf5files = glob.glob(os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"] + "*.hdf5"))
	hdf5files.sort() # sort the files by number
	print("Running Plot routines!")

	# get list of plot variables
	plot_variables = config["VariableData"].get("plot_variables", "").split()
	print("Number of plot variables: ", len(plot_variables))
	print("Variables to be plotted: ", config["VariableData"].get("plot_variables", ""))

	# setup the compute engine
	setup_engine()

	#iterate over the plot variables, and create plots for each hdf5
	for i in range(len(plot_variables)):
		plotvar = plot_variables[i]
		
		#get bounds for this plot variables
		config_setplotbounds_string = "set_plot_bounds_" + str(i+1)
		config_plotbounds_string = "plot_variables_range_" + str(i+1)
		print ("Plotting slices for {0}...".format(plotvar))
		setplotbounds = config["VariableData"].getint(config_setplotbounds_string, 0)
		plotbounds = tuple(np.float64(config["VariableData"].get(config_plotbounds_string, "0 0").split()))

		#create the slice plot
		make_slice_plots(plotvar, hdf5files, setplotbounds, plotbounds)

		#make a movie if asked
		if MultipleDatabase() & config["Output"].getint("make_movie", fallback = 0):
			print ("Making a movie...")
			cmd = 'ffmpeg -r ' + str(config["Output"].get("movie_framerate", fallback = "5")) + ' -s 1920x1080 -i ' + plotpath +"/"+ plotvar + '%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ' +moviepath +"/"+ plotvar+ '.mp4'
			system_status = os.system(cmd)
			
			#Check ffmpeg status to ensure command executed successfully
			if not system_status == 0:
				raise 	OSError("ffmpeg failed. Could not make the movie")

		print("I've finished!")
	
	os.remove("./visitlog.py")
	closedatabasesuccess = CloseDatabase()
	closeenginesuccess = CloseComputeEngine()

	if not closedatabasesuccess:
		raise IOError("Could not close the database")
	if not closeenginesuccess:
		raise IOError("Could not close the compute engine")

	exit()

if __visit_script_file__ == __visit_source_file__:
	main()
