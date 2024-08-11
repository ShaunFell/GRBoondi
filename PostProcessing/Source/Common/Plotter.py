## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

from .Utils import *
from Source.TwoD import *
from Source.ThreeD import *

@require_visit
def make_visit_plot(config, variableToPlot, hdf5files, plot_type = '2d', setplotbounds = False, plotbounds = [0,1]) :
	""" This function iterates over all hdf5 files and generates the plots according to the plot_type

	Args:
		config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters
		variableToPlot (str): name of the variable to be plotted
		hdf5files (list): list of paths to the hdf5 files
		plot_type (str, optional): type of plot to generate. Defaults to "2d". Options are "2d" and "3d"
		setplotbounds (bool, optional): flag to specify whether bounds for the plotting variable should be specified. Defaults to False.
		plotbounds (list, optional): list of two values that specify the minimum and maximum value of the plot variable. Defaults to [0,1].

	Raises:
		IOError: Database could not be opened
		RuntimeError: TimeSlider could not advance
		RuntimeError: Window could not be saved
		SystemError: Plot could not be cleared from memory
	"""
	import visit

	# open all the hdf5 files.
	# create a database object if there is more than one
	Open_Database(config)

	# suppress all messages except warnings and errors
	visit.SuppressMessages(config["Header"].getint("verbosity",0)) # https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/python_scripting/functions.html#suppressmessages		

	# setup verbosity printing	
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))
	verbPrint("Verbosity: ", verbPrint.verbosity)

	#Determine starting point, if overwrite deactivated
	for i in range(1, len(hdf5files)):
		firstfilename = str(variableToPlot) + ('%04d' % i)
		firstfilepath =  os.path.join(config["Output"]["output_plot_path"], firstfilename +"."+ config["Output"].get("fileform").lower())
		if not config["Output"].getboolean("overwrite_plots",0) and os.path.exists(firstfilepath):
			verbPrint("Plot already exists. Skipping...")
			timeslider_status = visit.TimeSliderNextState() # Advance to next state
			continue
		else:
			# file doesnt exist, so we should start the plotting here
			timeslider_status = visit.TimeSliderNextState()
			break	
	
		if not timeslider_status:
			raise RuntimeError("TimeSlider could not advance!")

	# create the plot object based on plot type
	plot_func_selector(plot_type)(config, variableToPlot, plotbounds, setplotbounds) #includes error handling

	
	# iterate over all hdf5 files, and create a new plot for each
	# then save to disk
	for i in range(1, len(hdf5files)):
		savename = str(variableToPlot) + ('%04d' % i)
		save_abs_path = os.path.join(config["Output"]["output_plot_path"], savename +"."+ config["Output"].get("fileform").lower())
		
		
		#if the plot already exists and overwrite disabled, skip
		if not config["Output"].getboolean("overwrite_plots",0) and os.path.exists(save_abs_path): 
			verbPrint("Plot already exists. Skipping...")
			visit.TimeSliderNextState() # Advance to next state
			continue

		verbPrint("Plotting file " + hdf5files[i])
		timeslider_status = visit.TimeSliderNextState() #advance to next state

		if not timeslider_status:
			raise RuntimeError("TimeSlider could not advance!")

		#save the window to file
		Save_Current_Window(config, variableToPlot)

	# clean up and close window
	plotdelete_status = visit.DeleteAllPlots()
	if not plotdelete_status:
		raise SystemError("Plot could not be deleted!")





