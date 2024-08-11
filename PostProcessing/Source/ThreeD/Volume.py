## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

from Source.Common.PlotSetup import *
from Source.Common.Utils import require_visit

@require_visit
def generate_3dslice_plot(config, variableToPlot, plotbounds = [0,1], setplotbounds = False):
	"""Generates a single volume plot for a given variable

	Args:
		config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters
		variableToPlot (str): name of the variable to be plotted
		setplotbounds (bool, optional): flag to specify whether bounds for the plotting variable should be specified. Defaults to False.
		plotbounds (list, optional): list of two values that specify the minimum and maximum value of the plot variable. Defaults to [0,1].

	Raises:
		RuntimeError: DrawPlots() failed. Aborting.
	"""

	# annotation settings, such as the information printed on the plot (user, database name, time, etc.)
	Annotation_Setup(config)

	## Add a volume plot
	Volume_Setup(config, variableToPlot)

	# plot all levels
	silr = visit.SILRestriction()
	silr.TurnOnAll()
	visit.SetPlotSILRestriction(silr,1)

	# Phew, now we can finally draw the plots!
	drawplot_status = visit.DrawPlots()
	if not drawplot_status:
		raise RuntimeError("DrawPlots() failed. Aborting.")

	## Set camera options
	View_3D_Setup(config)

	#Please save me! (to disk)
	Save_Current_Window(config, variableToPlot)
