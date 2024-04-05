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

def setup_slice_plot(variableToPlot, plotbounds, setplotbounds) :
	# annotation settings
	AnnotationAtts = AnnotationAttributes()
	AnnotationAtts.axes3D.visible = 1
	AnnotationAtts.userInfoFlag = 0
	AnnotationAtts.databaseInfoFlag = 1
	AnnotationAtts.timeInfoFlag = 1
	AnnotationAtts.legendInfoFlag = 1
	AnnotationAtts.axesArray.visible = 0
	AnnotationAtts.backgroundColor = (0, 0, 0, 255)
	AnnotationAtts.foregroundColor = (255, 255, 255, 255)
	AnnotationAtts.axes3D.bboxFlag = 1
	AnnotationAtts.axes3D.triadFlag = 0
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
	PseudocolorAtts.scaling = PseudocolorAtts.Linear # Linear, Log, Skew
	PseudocolorAtts.skewFactor = 1
     # limitsMode =  OriginalData, CurrentPlot
	PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData
	PseudocolorAtts.minFlag = setplotbounds
	PseudocolorAtts.min = plotbounds[0]
	PseudocolorAtts.maxFlag = setplotbounds
	PseudocolorAtts.max = setplotbounds[1]
	PseudocolorAtts.colorTableName = config["PlotConfig"].get("colortable", fallback = "viridis")
	PseudocolorAtts.invertColorTable = 0
     # opacityType = ColorTable, FullyOpaque, Constant, Ramp, VariableRange
	PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque 
	PseudocolorAtts.smoothingLevel = 0
     # centering = Natural, Nodal, Zonal
	PseudocolorAtts.centering = PseudocolorAtts.Nodal
	SetPlotOptions(PseudocolorAtts)
	
	# slice the pseudocolour plot
	AddOperator("Slice", 1)
	SliceAtts = SliceAttributes()
     # originType = Point, Intercept, Percent, Zone, Node
	SliceAtts.originType = SliceAtts.Point  
	originpoint = np.float64(config["PlotConfig"]["origin"].split())
	SliceAtts.originPoint = tuple(np.float64(config["PlotConfig"]["origin"].split()))
     # axisType = XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.axisType = SliceAtts.Arbitrary  
	SliceAtts.normal =tuple(np.float64(config["PlotConfig"]["normal_vec"].split()))
	SliceAtts.project2d = 1
	SliceAtts.flip = 0
	SetOperatorOptions(SliceAtts, 1)

	# plot all levels
	silr = SILRestriction()
	silr.TurnOnAll()
	SetPlotSILRestriction(silr ,1)
	
	# engage!
	DrawPlots()
	
	# Zoom in
	View2DAtts = View2DAttributes()
	View2DAtts.windowCoords = (
		config["PlotConfig"]["plot_u_min"], 
		config["PlotConfig"]["plot_u_max"],
		config["PlotConfig"]["plot_v_min"],
		config["PlotConfig"]["plot_v_max"]
	)
	View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
	View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
	View2DAtts.fullFrameAutoThreshold = 100
	View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
	View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
	View2DAtts.windowValid = 1
	SetView2D(View2DAtts)
	
	# Save Me
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = plotpath
	SaveWindowAtts.fileName = str(variableToPlot)
	SaveWindowAtts.family = 1
  # format = BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, 
  # RGB, STL, TIFF, ULTRA, VTK, PLY
	SaveWindowAtts.format = SaveWindowAtts.PNG
	SaveWindowAtts.width = 1024
	SaveWindowAtts.height = 1024
	SaveWindowAtts.quality = 80
	# resConstraint = NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.resConstraint = SaveWindowAtts.EqualWidthHeight 
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()


def PlotFiles():
	filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
	files = glob.glob(filename_prefix+ "*.3d.hdf5")
	plot_files = [x for x in files if config["Header"]["plot_header"] in x]
	return plot_files

def MultipleDatabase():
	
	if len(PlotFiles())>1:
		print("Opening multiple plot files")
		return True
	else:
		print("Opening single plot file")
		return False

def make_slice_plots(variableToPlot, hdf5files, setplotbounds, plotbounds) :
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
	hdf5files = glob.glob(config["Header"]["hdf5_path"] + config["Header"]["plot_header"] + "*.hdf5")
	hdf5files.sort()
	print("Running Plot routines!")

	plot_variables = config["VariableData"].get("plot_variables", "").split()
	print("Number of plot variables: ", len(plot_variables))
	print("Variables to be plotted: ", config["VariableData"].get("plot_variables", ""))
	for i in range(len(plot_variables)):
		plotvar = plot_variables[i]
		print ("Plotting slices for {0}...".format(plotvar))
		setplotbounds = config["VariableData"].get("set_plot_bounds_"+str(i), 0)
		plotbounds = np.float64(np.array(config["VariableData"].get("plot_bounds_"+str(i), "0 0").split()))
		make_slice_plots(plotvar, hdf5files, setplotbounds, plotbounds)

		if MultipleDatabase():
			print ("Making a movie...")
			os.system('ffmpeg -r 5 -f image2 -s 1920x1080 -i ' + plotpath + '%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ' +moviepath+ '.mp4')

		print("I've finished!")
	
	os.remove("./visitlog.py")
	exit()

if __visit_script_file__ == __visit_source_file__:
	main()
