#
## Fancier version of PlotSlice.py which makes a movie and uses
## time stepping to go through the files (thanks to Tiago for the tips)
## Run this using 'visit -nowin -cli -s PlotSliceRemix.py (path to hdf) (path to plot output)'
## to run in parallel do 'visit -np 32 -nowin -cli -s PlotSliceRemix.py'
## NB you will need ffmpeg installed to make the movie
#

import glob, os, sys

# --------------------------------------------------------------------------
## Inputs - you need to update these!
# --------------------------------------------------------------------------

# file details
plt_prefix = "GeneralizedProcap_"
path_to_hdf5_files = sys.argv[1]

# plot details
# select variable 
plot_variables = ["rho", "rhoE", "Ham", "Asquared"]


if len(sys.argv) >= 2:
	output_directory = sys.argv[2]
else:
	output_directory = os.getcwd()

if not os.path.exists(output_directory):
	os.mkdir(output_directory)
if not os.path.exists(output_directory+"../Movies"):
	os.mkdir(output_directory+"../Movies")

print("HDF5 file path: " + path_to_hdf5_files)
print("Output directory: " + output_directory)

# max and min values for colourbar
set_min_max = 0 # 1 for true, 0 for false
min_value = 0.0
max_value = 1.0
# slice origin and normal direction
origin_point_x = 32
origin_point_y = 32
origin_point_z = 0
normal_in_x = 0
normal_in_y = 0
normal_in_z = 1
# max and min coords in the sliced plane (e.g. x and y if normal to z)
# NB relative to origin as defined above
min_u = -20
max_u = 20
min_v = min_u
max_v = max_u

# --------------------------------------------------------------------------
# From here you only need to amend if you want to do something non standard
# --------------------------------------------------------------------------

def setup_slice_plot(variableToPlot) :
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
	PseudocolorAtts.minFlag = set_min_max
	PseudocolorAtts.min = min_value
	PseudocolorAtts.maxFlag = set_min_max
	PseudocolorAtts.max = max_value
	PseudocolorAtts.colorTableName = "viridis"
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
	SliceAtts.originPoint = (origin_point_x, origin_point_y, origin_point_z)
     # axisType = XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.axisType = SliceAtts.Arbitrary  
	SliceAtts.normal = (normal_in_x, normal_in_y, normal_in_z)
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
	View2DAtts.windowCoords = (min_u, max_u, min_v, max_v)
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
	SaveWindowAtts.outputDirectory = output_directory
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
	files = glob.glob(path_to_hdf5_files + "*.3d.hdf5")
	plot_files = [x for x in files if plt_prefix in x]
	return plot_files

def MultipleDatabase():
	
	if len(PlotFiles())>1:
		print("Opening multiple plot files")
		return True
	else:
		print("Opening single plot file")
		return False

def make_slice_plots(variableToPlot, hdf5files, hdf5files_base) :
	if MultipleDatabase():
		OpenDatabase(hdf5files_base + "*" + ".3d.hdf5 database", 0)
	else:
		OpenDatabase(PlotFiles()[0], 0)

	hdf5file = hdf5files[0]
	print("Setup and plot first slice from file : " + hdf5file)
	setup_slice_plot(variableToPlot)
	for i in range(1, len(hdf5files)):
		print("Plotting file " + hdf5files[i])
		TimeSliderNextState()
		SaveWindowAtts = SaveWindowAttributes()
		SaveWindowAtts.outputToCurrentDirectory = 0
		SaveWindowAtts.outputDirectory = output_directory
		SaveWindowAtts.fileName = str(variableToPlot) + ('%04d' % i)
		SaveWindowAtts.family = 0 # needs to be enforced again
		SetSaveWindowAttributes(SaveWindowAtts)	
		SaveWindow()
	DeleteAllPlots()



def main():
	hdf5files_base = path_to_hdf5_files + plt_prefix
	hdf5files = glob.glob(hdf5files_base + "*.hdf5")
	hdf5files.sort()
	print("Running Plot routines!")

	for plotvar in plot_variables:
		print ("Plotting slices for {0}...".format(plotvar))
		make_slice_plots(plotvar, hdf5files, hdf5files_base)

		if MultipleDatabase():
			print ("Making a movie...")
			os.system('ffmpeg -r 5 -f image2 -s 1920x1080 -i ' + output_directory+str(plotvar) + '%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p ' + output_directory+"../Movies/"+str(plotvar) + '.mp4')

		print("I've finished!")
	
	os.remove("./visitlog.py")
	exit()

if __visit_script_file__ == __visit_source_file__:
	main()
