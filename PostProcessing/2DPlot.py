#
## Run this using 'visit -nowin -cli -s 2DPlot.py paramfile'
## to run in parallel do 'visit -np 32 -nowin -cli -s 2DPlot.py paramfile'
## NOTE: you will need ffmpeg installed to make the movie
#

import glob, os, sys, configparser
import numpy as np

## Load config file
config = configparser.ConfigParser()
config.read(sys.argv[1])


# Extract config variables
plotpath = config["Output"]["output_plot_path"]
moviepath = config["Output"]["output_movie_path"]
overwrite_plots = config["Output"].getboolean("overwrite_plots",0)

#if the plot and movie paths dont exist, create them
if not os.path.exists(plotpath):
	os.mkdir(plotpath)
if not os.path.exists(moviepath):
	os.mkdir(moviepath)

print("HDF5 file path: " + config["Header"]["hdf5_path"])
print("Output Plots directory: " +plotpath)
print("Output Movie Directory: " +moviepath)

#set verbosity for printing
def verbPrint(*objects):
	if config["Header"].getint("verbosity",0):
		print(*objects)

verbPrint("Verbosity: ", config["Header"].getint("verbosity",0))


def setup_engine():
	"""
	Create the engine  for plotting
	"""
	useparallel = config["EngineConfig"].getboolean("use_parallel", 0)
	
	#if useparallel is disabled, default engine is launched
	if useparallel:
		host = config["EngineConfig"].get("host", "localhost")
		print("Running in parallel on host: ", host)
	
		# extract all the config options for the engine
		num_procs = config["EngineConfig"].get("number_processes", 1)
		num_nodes = config["EngineConfig"].get("number_nodes", 1)
		partition_name = config["EngineConfig"].get("partition")
		job_cmd = config["EngineConfig"].get("job_cmd", "srun")
		time_limit = config["EngineConfig"].get("time_limit", "01:00:00")
		job_name = config["EngineConfig"].get("job_name", "3DPlot")
		add_sub_args = config["EngineConfig"].get("additional_sub_args", "")

		arg = ("-l", job_cmd, "-n", job_name, "-p", partition_name, "-np", num_procs, "-nn", num_nodes, "-t", time_limit, "-la", add_sub_args)

		#Submit the job using visit job scheduler functionality
		OpenComputeEngine(host, arg)

# This is where the real magic happens. This function creates a single slice plot.
def setup_slice_plot(variableToPlot, plotbounds, setplotbounds) :
	"""
	Create all the data for a single slice plot
	"""

	"""
	For documentation on all the attributes below, see
		https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/python_scripting/functions.html
	"""


	# annotation settings, such as the information printed on the plot (user, database name, time, etc.)
	AnnotationAtts = AnnotationAttributes()
	axes3dvisible = config["AnnotationConfig"].getint('axes3Dvisible', 1)
	verbPrint("Axes3d visible: ", str(axes3dvisible))
	userinfoflag = config["AnnotationConfig"].getint('userInfoFlag', 0)
	verbPrint("userInfoFlag: ",str(userinfoflag))
	databaseinfoflag = config["AnnotationConfig"].getint('databaseInfoFlag', 1)
	verbPrint("databaseInfoFlag: " , str(databaseinfoflag))
	timeinfoflag = config["AnnotationConfig"].getint('timeInfoFlag', 1)
	verbPrint("timeInfoFlag: " , str(timeinfoflag))
	legendinfoflag = config["AnnotationConfig"].getint('legendInfoFlag', 1)
	verbPrint("legendInfoFlag: " , str(legendinfoflag))
	axesarrayvisible = config["AnnotationConfig"].getint('axesArrayvisible', 1)
	verbPrint('axesarrayvisible: ' , str(axesarrayvisible))
	backgroundcolor = tuple(map(int,config["AnnotationConfig"].get("backgroundColor", "0 0 255").split()))
	verbPrint("backgroundColor: " , backgroundcolor)
	foregroundcolor = tuple(map(int, config["AnnotationConfig"].get("foregroundColor", "255 255 255").split()))
	verbPrint("foregroundColor: " , foregroundcolor)
	bboxflag = config["AnnotationConfig"].getint('bboxFlag', 1)
	verbPrint("bboxFlag: " , str(bboxflag))
	triadflag = config["AnnotationConfig"].getint('triadFlag', 0)
	verbPrint("triadFlag: " , str(triadflag))

	#Set the annotation attributes
	AnnotationAtts.axes3D.visible = axes3dvisible
	AnnotationAtts.userInfoFlag = userinfoflag
	AnnotationAtts.databaseInfoFlag = databaseinfoflag
	AnnotationAtts.timeInfoFlag = timeinfoflag
	AnnotationAtts.legendInfoFlag = legendinfoflag
	AnnotationAtts.axesArray.visible = axesarrayvisible
	AnnotationAtts.backgroundColor = backgroundcolor
	AnnotationAtts.foregroundColor = foregroundcolor
	AnnotationAtts.axes3D.bboxFlag = bboxflag
	AnnotationAtts.axes3D.triadFlag = triadflag
	SetAnnotationAttributes(AnnotationAtts)
	
	# save settings
	InvertBackgroundColor()	
	
	# add pseudocolour plot
	AddPlot("Pseudocolor", variableToPlot, 1, 1)
	PseudocolorAtts = PseudocolorAttributes() #Get the pseudocolor plot attributes instance
	plotscaling = config["PlotConfig"].get("plotscaling", fallback = "Linear").lower() #Determine how we want the plot ticks to be scaled.
	if plotscaling == "linear":
		PseudocolorAtts.scaling = PseudocolorAtts.Linear
	elif plotscaling == "log":
		PseudocolorAtts.scaling = PseudocolorAtts.Log
	elif plotscaling == "skew":
		PseudocolorAtts.scaling = PseudocolorAtts.Skew
	else:
		print("Unknown plot scaling: " + plotscaling, "\n setting to default: linear")
		PseudocolorAtts.scaling = PseudocolorAtts.Linear

	PseudocolorAtts.skewFactor = 1
    # limitsMode =  OriginalData, CurrentPlot
	PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData
	verbPrint("Setting plot bounds: ", setplotbounds)
	verbPrint("Plot Bounds lower: ", plotbounds[0], " upper: ", plotbounds[1])
	PseudocolorAtts.minFlag = setplotbounds #Do we want a minimum value?
	PseudocolorAtts.min = plotbounds[0] #Set the minimum value
	PseudocolorAtts.maxFlag = setplotbounds #Do we want a maximum value?
	PseudocolorAtts.max = plotbounds[1] #Set the maximum value
	PseudocolorAtts.colorTableName = config["PlotConfig"].get("colortable", fallback = "viridis") #The type of color data.
	PseudocolorAtts.invertColorTable = config["PlotConfig"].getint("invert_color_table", fallback = 0) #Invert the color scale
    # opacityType = ColorTable, FullyOpaque, Constant, Ramp, VariableRange
	PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque 
	PseudocolorAtts.smoothingLevel = 0
    # centering = Natural, Nodal, Zonal
	centering = config["PlotConfig"].get("centering", fallback = "Natural")
	if centering == "Natural":
		PseudocolorAtts.centering = PseudocolorAtts.Natural
	elif centering == "Nodal":
		PseudocolorAtts.centering = PseudocolorAtts.Nodal
	elif centering == "Zonal":
		PseudocolorAtts.centering = PseudocolorAtts.Zonal

	#Set the options we just created above to the pseudocolor plot
	SetPlotOptions(PseudocolorAtts)
	

	# slice the pseudocolour plot
	AddOperator("Slice", 1)
	SliceAtts = SliceAttributes()

	origintype = config["SliceConfig"].get("origin_type", fallback = "Intercept") #Set how we want to slice the plot. Point, Intercept, Percent, Zone, Node
	verbPrint("Origin type: ", origintype)
	if origintype=="Point":
			SliceAtts.originType = SliceAtts.Point
	elif origintype	== "Intercept":
			SliceAtts.originType = SliceAtts.Intercept
	elif origintype == "Percent":
			SliceAtts.originType = SliceAtts.Percent
	elif origintype == "Zone":
			SliceAtts.originType = SliceAtts.Zone
	elif origintype == "Node":
			SliceAtts.originType = SliceAtts.Node

	originpoint = tuple(np.float64(config["SliceConfig"]["origin"].split()))
	verbPrint("Slice plane origin: ", originpoint)
	SliceAtts.originPoint = originpoint
     # axisType = XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	axistype = config["SliceConfig"].get("axistype", fallback = "ZAxis").lower()
	if axistype == "zaxis":
		SliceAtts.axisType = SliceAtts.ZAxis
	elif axistype == "yaxis":
		SliceAtts.axisType = SliceAtts.YAxis
	elif axistype == "xaxis":
		SliceAtts.axisType = SliceAtts.XAxis
	elif axistype == "arbitrary":
		SliceAtts.axisType = SliceAtts.Arbitrary
	elif axistype == "thetaphi":
		SliceAtts.axisType = SliceAtts.ThetaPhi

	#normal vector for the plane that slices the pseudocolor plot
	slicenormalvec = tuple(np.float64(config["SliceConfig"]["normal_vec"].split()))
	verbPrint("Normal vector: ", slicenormalvec)
	SliceAtts.normal = slicenormalvec #The slice plane's normal vector
	SliceAtts.project2d = 1 #Do we project the plot to 2D or keep the sliced surface in a 3D plot?
	SliceAtts.flip = 0
	SetOperatorOptions(SliceAtts, 1)


	#Should we plot the mesh as well?
	if config["MeshConfig"].getboolean("activatemesh", 0):
		AddPlot("Mesh", "Mesh",1,1)
		#Mesh attributes
		MeshAtts = MeshAttributes()
		MeshAtts.opacity = config["MeshConfig"].getfloat("meshopacity",1.0)
		meshcolor = tuple(np.int64(config["MeshConfig"]["meshcolor"].split()))
		verbPrint("Mesh color: ", meshcolor)
		MeshAtts.meshColor = meshcolor
		SetPlotOptions(MeshAtts)


	# plot all levels
	silr = SILRestriction()
	silr.TurnOnAll()
	SetPlotSILRestriction(silr ,1)
	
	# Ok, draw the plot now!
	DrawPlots()
	

	## Camera Settings

	# set the zoom
	View2DAtts = View2DAttributes()
	windowcoords = (
		config["ViewConfig"].getfloat("plot_u_min"),
		config["ViewConfig"].getfloat("plot_u_max"),
		config["ViewConfig"].getfloat("plot_v_min"),
		config["ViewConfig"].getfloat("plot_v_max") 
		)
	viewportcoords = tuple(np.float64(config["ViewConfig"]["viewportcoords"].split()))
	verbPrint("Window coordinates: ", windowcoords)
	verbPrint("View port coords: " , viewportcoords)
	
	View2DAtts.windowCoords = windowcoords
	View2DAtts.viewportCoords = viewportcoords
	View2DAtts.fullFrameActivationMode = View2DAtts.Off  # On, Off, Auto
	View2DAtts.fullFrameAutoThreshold = 100
	
	xscaling = config["PlotConfig"].get("xscaling", fallback = "Linear")
	yscaling = config["PlotConfig"].get("yscaling", fallback = "Linear")
	
	#set x and y scaling 
	if xscaling == "Linear":
		View2DAtts.xScale = View2DAtts.LINEAR 
	elif xscaling == "Log":
		View2DAtts.xScale = View2DAtts.LOG
	elif xscaling == "Skew":
		View2DAtts.xScale = View2DAtts.SKEW 
	if yscaling == "Linear":
		View2DAtts.yScale = View2DAtts.LINEAR
	elif yscaling == "Log":
		View2DAtts.yScale = View2DAtts.LOG
	elif yscaling == "Skew":
		View2DAtts.yScale = View2DAtts.SKEW

	View2DAtts.windowValid = 1
	SetView2D(View2DAtts)
	
	# Please save me! (to disk)
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = plotpath
	SaveWindowAtts.fileName = str(variableToPlot)
	SaveWindowAtts.family = 1
  
	# Set the file format
	outformat = config["Output"].get("fileform", fallback = "PNG").upper()
	if outformat == "png":
		SaveWindowAtts.format = SaveWindowAtts.PNG
	elif outformat == "JPG":
		SaveWindowAtts.format = SaveWindowAtts.JPG
	elif outformat == "TIF":
		SaveWindowAtts.format = SaveWindowAtts.TIFF
	elif outformat == "BMP":
		SaveWindowAtts.format = SaveWindowAtts.BMP
	elif outformat == "PS":
		SaveWindowAtts.format = SaveWindowAtts.POSTSCRIPT
	elif outformat == "EPS":
		SaveWindowAtts.format = SaveWindowAtts.POSTSCRIPT
	elif outformat == "PPM":
		SaveWindowAtts.format = SaveWindowAtts.PPM
	elif outformat == "RGB":
		SaveWindowAtts.format = SaveWindowAtts.RGB
	elif outformat == "STL":
		SaveWindowAtts.format = SaveWindowAtts.STL
	elif outformat == "EXR":
		SaveWindowAtts.format = SaveWindowAtts.EXR
	elif outformat == "ULTRA":
		SaveWindowAtts.format = SaveWindowAtts.ULTRA
	elif outformat == "VTK":
		SaveWindowAtts.format = SaveWindowAtts.VTK
	elif outformat == "PLY":
		SaveWindowAtts.format = SaveWindowAtts.PLY


	# Set the resolution
	SaveWindowAtts.width = np.float64(config["Output"].get("width", fallback = 1024))
	SaveWindowAtts.height = np.float64(config["Output"].get("height", fallback = 1024))
	SaveWindowAtts.quality = np.float64(config["Output"].get("quality", fallback = 80))
	# resConstraint = NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint 
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()


def PlotFiles():
	"""
	Find all the plot files and return them as list of absolute path strings
	"""
	filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
	files = glob.glob(filename_prefix+ "*.3d.hdf5")
	plot_files = [x for x in files if config["Header"]["plot_header"] in x]
	if len(plot_files) == 0:
		raise SystemError("No Plot Files Found!")
	
	return plot_files

def MultipleDatabase():
	""" 
	Flag that tells us if there's multiple files
	"""
	
	if len(PlotFiles())>1:
		print("Opening multiple plot files")
		return True
	else:
		print("Opening single plot file")
		return False

def make_slice_plots(variableToPlot, hdf5files, setplotbounds, plotbounds) :
	""" 
	Do that actual plotting, iterating over all the files for a given variable
	"""

	# open all the hdf5 files.
	# create a database object if there is more than one
	if MultipleDatabase():
		filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
		OpenDatabase(filename_prefix + "*" + ".3d.hdf5 database", 0)
	else:
		OpenDatabase(PlotFiles()[0], 0)

	#Determine starting point, if overwrite deactivated
	for i in range(1, len(hdf5files)):
		firstfilename = str(variableToPlot) + ('%04d' % i)
		firstfilepath =  os.path.join(plotpath, firstfilename +"."+ config["Output"].get("fileform").lower())
		if not overwrite_plots and os.path.exists(firstfilepath):
			verbPrint("Plot already exists. Skipping...")
			TimeSliderNextState() # Advance to next state
			continue
		else:
			# file doesnt exist, so we should start the plotting here
			TimeSliderNextState()
			break

	# create the plot
	setup_slice_plot(variableToPlot, plotbounds, setplotbounds)
	
	# iterate over all hdf5 files, and create a new plot for each
	# then save to disk
	for i in range(1, len(hdf5files)):
		savename = str(variableToPlot) + ('%04d' % i)
		save_abs_path = os.path.join(plotpath, savename +"."+ config["Output"].get("fileform").lower())
		
		
		#if the plot already exists and overwrite disabled, skip
		if not overwrite_plots and os.path.exists(save_abs_path): 
			verbPrint("Plot already exists. Skipping...")
			TimeSliderNextState() # Advance to next state
			continue

		print("Plotting file " + hdf5files[i])
		TimeSliderNextState() #advance to next state

		#save the window to file
		SaveWindowAtts = SaveWindowAttributes()
		SaveWindowAtts.outputToCurrentDirectory = 0
		SaveWindowAtts.outputDirectory = config["Output"]["output_plot_path"]
		SaveWindowAtts.fileName = str(variableToPlot) + ('%04d' % i)
		SaveWindowAtts.family = 0 # needs to be enforced again
		SetSaveWindowAttributes(SaveWindowAtts)	
		SaveWindow()

	# clean up and close window
	DeleteAllPlots()



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
			cmd = 'ffmpeg -r ' + str(config["Output"].get("movie_framerate", fallback = "5")) + '-s 1920x1080 -i ' + plotpath +"/"+ plotvar + '%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ' +moviepath +"/"+ plotvar+ '.mp4'
			os.system(cmd)

		print("I've finished!")
	
	os.remove("./visitlog.py")
	exit()

if __visit_script_file__ == __visit_source_file__:
	main()
