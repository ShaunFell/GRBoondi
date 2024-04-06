#
## Run this using 'visit -nowin -cli -s PlotSliceRemix.py (path to hdf) (param file) (optional output directory)'
## to run in parallel do 'visit -np 32 -nowin -cli -s PlotSliceRemix.py'
## NOTE: you will need ffmpeg installed to make the movie
#

import glob, os, sys, configparser
import numpy as np

## Load config file
config = configparser.ConfigParser()
config.read(sys.argv[1])



plotpath = config["Output"]["output_plot_path"]
moviepath = config["Output"]["output_movie_path"]

if not os.path.exists(plotpath):
	os.mkdir(plotpath)
if not os.path.exists(moviepath):
	os.mkdir(moviepath)

print("HDF5 file path: " + config["Header"]["hdf5_path"])
print("Output Plots directory: " +plotpath)
print("Output Movie Directory: " +moviepath)

def verbPrint(*objects):
	if config["Header"].getint("verbosity",0):
		print(*objects)

verbPrint("Verbosity: ", config["Header"].getint("verbosity",0))

def setup_slice_plot(variableToPlot, plotbounds, setplotbounds) :
	# annotation settings
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
	saveAtts = SaveWindowAttributes()
	saveAtts.format = saveAtts.PNG
	saveAtts.resConstraint = saveAtts.NoConstraint
	saveAtts.width = 1280
	saveAtts.height = 620
	SetSaveWindowAttributes(saveAtts) 
	
	
	# add pseudocolour plot
	AddPlot("Pseudocolor", variableToPlot, 1, 1)
	PseudocolorAtts = PseudocolorAttributes()
	plotscaling = config["PlotConfig"].get("plotscaling", fallback = "Linear").lower()
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
	PseudocolorAtts.minFlag = setplotbounds
	PseudocolorAtts.min = plotbounds[0]
	PseudocolorAtts.maxFlag = setplotbounds
	PseudocolorAtts.max = plotbounds[1]
	PseudocolorAtts.colorTableName = config["PlotConfig"].get("colortable", fallback = "viridis") #Color table name
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

	SetPlotOptions(PseudocolorAtts)
	
	# slice the pseudocolour plot
	AddOperator("Slice", 1)
	SliceAtts = SliceAttributes()

	origintype = config["SliceConfig"].get("origin_type", fallback = "Intercept")
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

	#normal vector for slice surface
	slicenormalvec = tuple(np.float64(config["SliceConfig"]["normal_vec"].split()))
	verbPrint("Normal vector: ", slicenormalvec)
	SliceAtts.normal = slicenormalvec
	SliceAtts.project2d = 1
	SliceAtts.flip = 0
	SetOperatorOptions(SliceAtts, 1)

	# plot all levels
	silr = SILRestriction()
	silr.TurnOnAll()
	SetPlotSILRestriction(silr ,1)
	
	# Ok, draw the plot now!
	DrawPlots()
	
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
  
	outformat = config["PlotConfig"].get("format", fallback = "PNG").upper()
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
	return plot_files

def MultipleDatabase():
	""" Flag that tells us if theres multiple files
	"""
	
	if len(PlotFiles())>1:
		print("Opening multiple plot files")
		return True
	else:
		print("Opening single plot file")
		return False

def make_slice_plots(variableToPlot, hdf5files, setplotbounds, plotbounds) :
	""" Do that actual plotting, iterating over all the files for a given variable
	"""
	if MultipleDatabase():
		filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
		OpenDatabase(filename_prefix + "*" + ".3d.hdf5 database", 0)
	else:
		OpenDatabase(PlotFiles()[0], 0)

	hdf5file = hdf5files[0]
	print("Setup and plot first slice from file : " + hdf5file)
	setup_slice_plot(variableToPlot, plotbounds, setplotbounds)
	for i in range(1, len(hdf5files)):
		print("Plotting file " + hdf5files[i])
		TimeSliderNextState()
		SaveWindowAtts = SaveWindowAttributes()
		SaveWindowAtts.outputToCurrentDirectory = 0
		SaveWindowAtts.outputDirectory = config["Output"]["output_plot_path"]
		SaveWindowAtts.fileName = str(variableToPlot) + ('%04d' % i)
		SaveWindowAtts.family = 0 # needs to be enforced again
		SetSaveWindowAttributes(SaveWindowAtts)	
		SaveWindow()
	DeleteAllPlots()



def main():

	hdf5files = glob.glob(os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"] + "*.hdf5"))
	hdf5files.sort()
	print("Running Plot routines!")

	plot_variables = config["VariableData"].get("plot_variables", "").split()
	print("Number of plot variables: ", len(plot_variables))
	print("Variables to be plotted: ", config["VariableData"].get("plot_variables", ""))
	for i in range(len(plot_variables)):
		plotvar = plot_variables[i]
		config_setplotbounds_string = "set_plot_bounds_" + str(i+1)
		config_plotbounds_string = "plot_variables_range_" + str(i+1)
		print ("Plotting slices for {0}...".format(plotvar))
		setplotbounds = config["VariableData"].getint(config_setplotbounds_string, 0)
		plotbounds = tuple(np.float64(config["VariableData"].get(config_plotbounds_string, "0 0").split()))
		make_slice_plots(plotvar, hdf5files, setplotbounds, plotbounds)

		if MultipleDatabase() & config["Output"].getint("make_movie", fallback = 0):
			print ("Making a movie...")
			os.system('ffmpeg -r 5 -f image2 -s 1920x1080 -i ' + plotpath + '%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ' +moviepath+ '.mp4')

		print("I've finished!")
	
	os.remove("./visitlog.py")
	exit()

if __visit_script_file__ == __visit_source_file__:
	main()
