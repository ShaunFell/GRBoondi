## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

from Utils import *
from 2D import *
from 3D import *


def make_slice_plots(config, variableToPlot, hdf5files, setplotbounds = False, plotbounds = [0,1], slice_type = "2d") :
    """ This function iterates over all hdf5 files and generates the plots. 

    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters
        variableToPlot (str): name of the variable to be plotted
        hdf5files (list): list of paths to the hdf5 files
        setplotbounds (bool, optional): flag to specify whether bounds for the plotting variable should be specified. Defaults to False.
        plotbounds (list, optional): list of two values that specify the minimum and maximum value of the plot variable. Defaults to [0,1].
        slice_type (str, optional): type of plot to generate. Defaults to "2d". Options are "2d" and "3d"

    Raises:
        IOError: Database could not be opened
        RuntimeError: TimeSlider could not advance
        RuntimeError: Window could not be saved
        SystemError: Plot could not be cleared from memory
    """

	# open all the hdf5 files.
	# create a database object if there is more than one
	if MultipleDatabase():
		filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
		database_status = OpenDatabase(filename_prefix + "*" + ".3d.hdf5 database", 0)
	else:
		database_status = OpenDatabase(PlotFiles()[0], 0)
	
	if not database_status:
		raise IOError("Could not open database!")

	#Determine starting point, if overwrite deactivated
	for i in range(1, len(hdf5files)):
		firstfilename = str(variableToPlot) + ('%04d' % i)
		firstfilepath =  os.path.join(plotpath, firstfilename +"."+ config["Output"].get("fileform").lower())
		if not config["Output"].getboolean("overwrite_plots",0) and os.path.exists(firstfilepath):
			verbPrint("Plot already exists. Skipping...")
			timeslider_status = TimeSliderNextState() # Advance to next state
			continue
		else:
			# file doesnt exist, so we should start the plotting here
			timeslider_status = TimeSliderNextState()
			break
	
		if not timeslider_status:
			raise RuntimeError("TimeSlider could not advance!")

	# create the plot
	setup_slice_plot(variableToPlot, plotbounds, setplotbounds)
	
	# iterate over all hdf5 files, and create a new plot for each
	# then save to disk
	for i in range(1, len(hdf5files)):
		savename = str(variableToPlot) + ('%04d' % i)
		save_abs_path = os.path.join(plotpath, savename +"."+ config["Output"].get("fileform").lower())
		
		
		#if the plot already exists and overwrite disabled, skip
		if not config["Output"].getboolean("overwrite_plots",0) and os.path.exists(save_abs_path): 
			verbPrint("Plot already exists. Skipping...")
			TimeSliderNextState() # Advance to next state
			continue

		print("Plotting file " + hdf5files[i])
		timeslider_status = TimeSliderNextState() #advance to next state

		if not timeslider_status:
			raise RuntimeError("TimeSlider could not advance!")

		#save the window to file
		SaveWindowAtts = SaveWindowAttributes()
		SaveWindowAtts.outputToCurrentDirectory = 0
		SaveWindowAtts.outputDirectory = config["Output"]["output_plot_path"]
		SaveWindowAtts.fileName = str(variableToPlot) + ('%04d' % i)
		SaveWindowAtts.family = 0 # needs to be enforced again
		SetSaveWindowAttributes(SaveWindowAtts)	
		windowsave_status = SaveWindow()
		if not windowsave_status:
			raise RuntimeError("Window could not be saved!")

	# clean up and close window
	plotdelete_status = DeleteAllPlots()
	if not plotdelete_status:
		raise SystemError("Plot could not be deleted!")





