## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

from Source.Common.PlotSetup import *
from Source.Common.Utils import require_visit

@require_visit
def generate_2dslice_plot(config, variableToPlot, plotbounds = [0,1], setplotbounds = False) :
	"""Generates a single slice plot for a given variable

	Args:
		config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters
		variableToPlot (str): name of the variable to be plotted
		setplotbounds (bool, optional): flag to specify whether bounds for the plotting variable should be specified. Defaults to False.
		plotbounds (list, optional): list of two values that specify the minimum and maximum value of the plot variable. Defaults to [0,1].

	Raises:
		RuntimeError: DrawPlots() failed. Aborting.
	"""

	## For documentation on all the attributes below, see
	## https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/python_scripting/functions.html

	# Set the annotation settings
	Annotation_Setup(config)
	
	# save settings
	visit.InvertBackgroundColor()	
	
	#Create a pseudocolor plot
	Pseudocolor_Setup(config, variableToPlot, setplotbounds, plotbounds)
	
	#Slice the pseudocolor plot
	SliceAttribute_Setup(config)
	
	#Should we plot the mesh as well?
	if config["MeshConfig"].getboolean("activatemesh", 0):
		Mesh_Setup(config)

	# plot all levels
	silr = visit.SILRestriction()
	silr.TurnOnAll()
	visit.SetPlotSILRestriction(silr ,1)
	
	# Ok, draw the plot now!
	drawplot_status = visit.DrawPlots()
	if not drawplot_status:
		raise RuntimeError("Could not draw plots!")

	# set the zoom
	View_2D_Setup(config)
	
	# Please save me! (to disk)
	Save_Current_Window(config, variableToPlot)
