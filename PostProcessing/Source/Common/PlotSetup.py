## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory


## For documentation on all the attributes below, see
## https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/python_scripting/functions.html
	
from Source.Common.Utils import *
import numpy as np

@require_visit
def Annotation_Setup(config):
	"""Set the annotation settings. 

	Args:
		config (_type_): _description_
	"""

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	# annotation settings, such as the information printed on the plot (user, database name, time, etc.)
	AnnotationAtts = visit.AnnotationAttributes()

	## Grab annotation settings from config file, with defaults
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
	visit.SetAnnotationAttributes(AnnotationAtts)

@require_visit
def Pseudocolor_Setup(config, variableToPlot, setplotbounds, plotbounds):

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	# add pseudocolour plot
	visit.AddPlot("Pseudocolor", variableToPlot, 1, 1)
	PseudocolorAtts = visit.PseudocolorAttributes() #Get the pseudocolor plot attributes instance
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
	visit.SetPlotOptions(PseudocolorAtts)

@require_visit
def Volume_Setup(config, variableToPlot):

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	visit.AddPlot("Volume", variableToPlot,1,1)
	VolumeAtts = visit.VolumeAttributes()
	plotscaling = config["PlotConfig"].get("plotscaling", fallback = "Linear").title()
	verbPrint("plot scaling: ",  plotscaling)

	#Set the attributes for the volume plot
	VolumeAtts.scaling = getattr(VolumeAtts, plotscaling) #Get the plot scaling attribute
	VolumeAtts.lightingFlag = config["PlotConfig"].getboolean('lightingFlag', 1)
	VolumeAtts.legendFlag =   config["PlotConfig"].getboolean('legendFlag', 1)
	VolumeAtts.opacityAttenuation = config["PlotConfig"].getint('opacityAttenuation', 1)
	opacitymode = config["PlotConfig"].get("opacityMode", fallback = "freeform").title()
	verbPrint("opacityMode: ",  opacitymode)
	VolumeAtts.opacityMode  = getattr(VolumeAtts, opacitymode+"Mode") #get the opacity mode attribute
	
	# add opacity scaling for variable data
	## set domain of width 256 and height 256-1
	domain = np.linspace(0,255,256) 
	#create the opacity using the piecewise function
	low = 255*config["PlotConfig"].getfloat("ramp_min", 0)
	high = 255*config["PlotConfig"].getfloat("ramp_max", 1)
	# add opacity ramp to the volume plot
	opacityramp =  tuple(np.piecewise(domain, [domain<low, ((domain>=low) & (domain<=high)), domain>high], [0, lambda x: 255*(x - low)/(high-low), 255]))
	verbPrint("opacity ramp: ",  opacityramp)
	VolumeAtts.freeformOpacity = opacityramp #set the tuple as the opacity ramp

	#Determine the max and min values of the variable to use for coloring and opacity
	VolumeAtts.useColorVarMin = config["PlotConfig"].getboolean('useColorVarMin', 0)
	VolumeAtts.colorVarMin = config["PlotConfig"].getfloat('colorVarMin', 0)
	VolumeAtts.useColorVarMax = config["PlotConfig"].getboolean('useColorVarMax', 0)
	VolumeAtts.colorVarMax = config["PlotConfig"].getfloat('colorVarMax', 0)
	VolumeAtts.useOpacityVarMin = config["PlotConfig"].getboolean('useOpacityVarMin', 0)
	VolumeAtts.opacityVarMin = config["PlotConfig"].getfloat('opacityVarMin', 0)
	VolumeAtts.useOpacityVarMax = config["PlotConfig"].getboolean('useOpacityVarMax', 0)
	VolumeAtts.opacityVarMax = config["PlotConfig"].getfloat('opacityVarMax', 0)

	#Specify the type of rendering engine to use
	rendertype = config["PlotConfig"].get("rendererType", fallback = "default").title()
	if rendertype == "Raycasting": rendertype = "RayCasting" 
	if rendertype == "Raycastingintegration": rendertype = "RayCastingIntegration"
	if rendertype == "Raycastingslivr": rendertype = "RayCastingSLIVR" 
	if rendertype == "Raycastingospray": rendertype = "RayCastingOSPRay" 
	verbPrint("Renderer type: ", rendertype)
	VolumeAtts.rendererType = getattr(VolumeAtts, rendertype) #get the rendering type attribute

	#Determine the sampling rate of the data to use for the plot
	sampling = config["PlotConfig"].get("sampling", fallback = "rasterization").title()
	verbPrint("Sampling: ", sampling)
	VolumeAtts.sampling = getattr(VolumeAtts, sampling) #get the sampling attribute
	lowgradientlightingreduc = config["PlotConfig"].get("lowGradientLightingReduction", fallback = "Lower").title()
	verbPrint("Low gradient lighting reduction: ", lowgradientlightingreduc)
	VolumeAtts.lowGradientLightingReduction = getattr(VolumeAtts, lowgradientlightingreduc) #get the low gradient lighting reduction attribute
	VolumeAtts.samplesPerRay = config["PlotConfig"].getint('samplesPerRay', 1)
	#set the above attributes to the plot options
	visit.SetPlotOptions(VolumeAtts)

@require_visit
def SliceAttribute_Setup(config):

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	# slice the pseudocolour plot
	visit.AddOperator("Slice", 1)
	SliceAtts = visit.SliceAttributes()

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
	visit.SetOperatorOptions(SliceAtts, 1)

@require_visit
def Mesh_Setup(config):

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	visit.AddPlot("Mesh", "Mesh",1,1)
	
	#Mesh attributes
	MeshAtts = visit.MeshAttributes()
	MeshAtts.opacity = config["MeshConfig"].getfloat("meshopacity",1.0)
	meshcolor = tuple(np.int64(config["MeshConfig"]["meshcolor"].split()))
	verbPrint("Mesh color: ", meshcolor)
	MeshAtts.meshColor = meshcolor
	visit.SetPlotOptions(MeshAtts)

@require_visit
def View_2D_Setup(config):

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	View2DAtts = visit.View2DAttributes()
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
	visit.SetView2D(View2DAtts)	

@require_visit
def View_3D_Setup(config):

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	View3DAtts = visit.View3DAttributes()
	viewnormal = tuple(np.float64(config["ViewConfig"].get("viewNormal", "0 0 1").split()))
	focus = tuple(np.float64(config["ViewConfig"].get("focus", "0 0 0").split()))
	viewup = tuple(np.float64(config["ViewConfig"].get("viewUp", "0 1 0").split()))
	View3DAtts.viewNormal = viewnormal
	View3DAtts.focus = focus
	View3DAtts.viewUp = viewup
	View3DAtts.viewAngle = config["ViewConfig"].getint("viewAngle", 30)
	View3DAtts.parallelScale = config["ViewConfig"].getfloat("parallelScale", 1.0)
	View3DAtts.nearPlane = config["ViewConfig"].getfloat("nearPlane", -0.1)
	View3DAtts.farPlane = config["ViewConfig"].getfloat("farPlane", 0.1)
	View3DAtts.imagePan = tuple(np.float64(config["ViewConfig"].get("imagePan", "0 0").split()))
	View3DAtts.imageZoom = config["ViewConfig"].getfloat("imageZoom", 1.0)
	View3DAtts.perspective = config["ViewConfig"].getint("perspective", 1)
	View3DAtts.eyeAngle = config["ViewConfig"].getint("eyeAngle", 0)
	View3DAtts.centerOfRotationSet = config["ViewConfig"].getint("centerOfRotationSet", 0)
	View3DAtts.centerOfRotation = tuple(np.float64(config["ViewConfig"].get("centerOfRotation", "0 0 0").split()))
	View3DAtts.axis3DScaleFlag = config["ViewConfig"].getboolean("axis3DScaleFlag", 0)
	View3DAtts.axis3DScales = tuple(np.float64(config["ViewConfig"].get("axis3DScales", "1 1 1").split()))
	View3DAtts.shear = tuple(np.float64(config["ViewConfig"].get("shear", "0 0 0").split()))
	View3DAtts.windowValid = config["ViewConfig"].getboolean("windowValid", 1)
	visit.SetView3D(View3DAtts)

@require_visit
def Save_Current_Window(config, variableToPlot):

	#initialize verbosity printing
	verbPrint = VerbosityPrint(config["Header"].getint("verbosity",0))

	SaveWindowAtts = visit.SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = config["Output"]["output_plot_path"]
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


	# Set the resolution, defaults to 1024x1024
	SaveWindowAtts.width = np.float64(config["Output"].get("width", fallback = 1024))
	SaveWindowAtts.height = np.float64(config["Output"].get("height", fallback = 1024))
	SaveWindowAtts.quality = np.float64(config["Output"].get("quality", fallback = 80))
	
	# No constraints on resolution, only user-specified
	# resConstraint = NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint 
	
	#set the window attributes just defined
	visit.SetSaveWindowAttributes(SaveWindowAtts)

	#finally, save the window
	visit.SaveWindow()
