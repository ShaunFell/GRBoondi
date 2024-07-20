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

		#Parse config options for engine
		num_procs = config["EngineConfig"].get("number_processes", 1)
		num_nodes = config["EngineConfig"].get("number_nodes", 1)
		partition_name = config["EngineConfig"].get("partition")
		job_cmd = config["EngineConfig"].get("job_cmd", "srun")
		time_limit = config["EngineConfig"].get("time_limit", "01:00:00")
		job_name = config["EngineConfig"].get("job_name", "3DPlot")
		add_sub_args = config["EngineConfig"].get("additional_launch_args", "")
		use_gpus = config["EngineConfig"].getboolean("use_gpus", 0)
		ngpus_per_node = config["EngineConfig"].get("ngpus_per_node", "1")

		#create argument tuple
		arg = ("-l", job_cmd, "-n", job_name, "-p", partition_name, "-np", num_procs, "-nn", num_nodes, "-t", time_limit, "-la", add_sub_args)

		#if gpus requested, add appropriate arguments to list
		if use_gpus:
			slurm_gpu_submission = "--gres=gpu:{0}".format(ngpus_per_node)
			add_sub_args += " {0}".format(slurm_gpu_submission)
			arg = arg[:-2] + ("-la", add_sub_args,)
			arg = arg + ("-hw_accel",)
			arg = arg + ("-n-gpus-per-node", ngpus_per_node,)
		

		#Submit the job
		openengine_status = OpenComputeEngine(host, arg)
		if not openengine_status:
			print("Job submission failed. Exiting...")
			sys.exit(1)



def setup_slice_plot(variableToPlot, plotbounds, setplotbounds):
	"""
	This is where the real magic happens. Create all the data for the plot
	"""

	"""
	For documentation on all the attributes below, see
		https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/python_scripting/functions.html
	"""

	# annotation settings, such as the information printed on the plot (user, database name, time, etc.)
	AnnotationAtts = AnnotationAttributes()

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

	## Set the annotation attributes
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



	## Add a volume plot

	AddPlot("Volume", variableToPlot,1,1)
	VolumeAtts = VolumeAttributes()
	plotscaling = config["VolumeConfig"].get("plotscaling", fallback = "Linear").title()
	verbPrint("plot scaling: ",  plotscaling)

	#Set the attributes for the volume plot
	VolumeAtts.scaling = getattr(VolumeAtts, plotscaling) #Get the plot scaling attribute
	VolumeAtts.lightingFlag = config["VolumeConfig"].getboolean('lightingFlag', 1)
	VolumeAtts.legendFlag =   config["VolumeConfig"].getboolean('legendFlag', 1)
	VolumeAtts.opacityAttenuation = config["VolumeConfig"].getint('opacityAttenuation', 1)
	opacitymode = config["VolumeConfig"].get("opacityMode", fallback = "freeform").title()
	verbPrint("opacityMode: ",  opacitymode)
	VolumeAtts.opacityMode  = getattr(VolumeAtts, opacitymode+"Mode") #get the opacity mode attribute
	
	# add opacity scaling for variable data
	## set domain of width 256 and height 256-1
	domain = np.linspace(0,255,256) 
	#create the opacity using the piecewise function
	low = 255*config["VolumeConfig"].getfloat("ramp_min", 0)
	high = 255*config["VolumeConfig"].getfloat("ramp_max", 1)
	# add opacity ramp to the volume plot
	opacityramp =  tuple(np.piecewise(domain, [domain<low, ((domain>=low) & (domain<=high)), domain>high], [0, lambda x: 255*(x - low)/(high-low), 255]))
	verbPrint("opacity ramp: ",  opacityramp)
	VolumeAtts.freeformOpacity = opacityramp #set the tuple as the opacity ramp

	#Determine the max and min values of the variable to use for coloring and opacity
	VolumeAtts.useColorVarMin = config["VolumeConfig"].getboolean('useColorVarMin', 0)
	VolumeAtts.colorVarMin = config["VolumeConfig"].getfloat('colorVarMin', 0)
	VolumeAtts.useColorVarMax = config["VolumeConfig"].getboolean('useColorVarMax', 0)
	VolumeAtts.colorVarMax = config["VolumeConfig"].getfloat('colorVarMax', 0)
	VolumeAtts.useOpacityVarMin = config["VolumeConfig"].getboolean('useOpacityVarMin', 0)
	VolumeAtts.opacityVarMin = config["VolumeConfig"].getfloat('opacityVarMin', 0)
	VolumeAtts.useOpacityVarMax = config["VolumeConfig"].getboolean('useOpacityVarMax', 0)
	VolumeAtts.opacityVarMax = config["VolumeConfig"].getfloat('opacityVarMax', 0)

	#Specify the type of rendering engine to use
	rendertype = config["VolumeConfig"].get("rendererType", fallback = "default").title()
	if rendertype == "Raycasting": rendertype = "RayCasting" 
	if rendertype == "Raycastingintegration": rendertype = "RayCastingIntegration"
	if rendertype == "Raycastingslivr": rendertype = "RayCastingSLIVR" 
	if rendertype == "Raycastingospray": rendertype = "RayCastingOSPRay" 
	verbPrint("Renderer type: ", rendertype)
	VolumeAtts.rendererType = getattr(VolumeAtts, rendertype) #get the rendering type attribute

	#Determine the sampling rate of the data to use for the plot
	sampling = config["VolumeConfig"].get("sampling", fallback = "rasterization").title()
	verbPrint("Sampling: ", sampling)
	VolumeAtts.sampling = getattr(VolumeAtts, sampling) #get the sampling attribute
	lowgradientlightingreduc = config["VolumeConfig"].get("lowGradientLightingReduction", fallback = "Lower").title()
	verbPrint("Low gradient lighting reduction: ", lowgradientlightingreduc)
	VolumeAtts.lowGradientLightingReduction = getattr(VolumeAtts, lowgradientlightingreduc) #get the low gradient lighting reduction attribute
	VolumeAtts.samplesPerRay = config["VolumeConfig"].getint('samplesPerRay', 1)
	#set the above attributes to the plot options
	SetPlotOptions(VolumeAtts)


	# plot all levels
	silr = SILRestriction()
	silr.TurnOnAll()
	SetPlotSILRestriction(silr,1)

	# Phew, now we can finally draw the plots!
	drawplot_status = DrawPlots()
	if not drawplot_status:
		raise Exception("DrawPlots() failed. Aborting.")

	## Set camera options

	# set the zoom
	View3DAtts = View3DAttributes()
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
	SetView3D(View3DAtts)



	#Please save me! (to disk)
	InvertBackgroundColor()
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = plotpath
	SaveWindowAtts.fileName = str(variableToPlot)
	SaveWindowAtts.family = 1

	#Set the output format
	outformat = config["Output"].get("fileform", fallback = "PNG").upper()
	verbPrint("File format: ", outformat)
	SaveWindowAtts.format = getattr(SaveWindowAtts, outformat)

	#Set the resolution
	SaveWindowAtts.width = np.float64(config["Output"].get("width", fallback = 1024))
	SaveWindowAtts.height = np.float64(config["Output"].get("height", fallback = 1024))
	SaveWindowAtts.quality = np.float64(config["Output"].get("quality", fallback = 80))
	# resConstraint = NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint 

	#save the window attributes
	SetSaveWindowAttributes(SaveWindowAtts)

	#execute the window save
	SaveWindow()




def PlotFiles():
	"""
	Find all the plot files and return them as list of absolute path strings
	"""
	filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
	files = glob.glob(filename_prefix+ "*.3d.hdf5")
	plot_files = [x for x in files if config["Header"]["plot_header"] in x]
	if len(plot_files) == 0:
		raise FileNotFoundError("No Plot Files Found!")
	return plot_files


def MultipleDatabase():
	""" 
	Flag that tells us if theres multiple files
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
		if config["Header"]["use_plot_range"]:
			timeindex = config["Header"].getint("plot_range", 0)
		else:
			timeindex = 0
		filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
		database_status = OpenDatabase(filename_prefix + "*" + ".3d.hdf5 database", timeindex)
	else:
		database_status = OpenDatabase(PlotFiles()[0], 0)
	
	#check database successfully opened
	if not database_status:
		raise SystemError("Database could not be opened!")

	#Determine starting point, if overwrite deactivated
	for i in range(1, len(hdf5files)):
		firstfilename = str(variableToPlot) + ('%04d' % i)
		firstfilepath =  os.path.join(plotpath, firstfilename +"."+ config["Output"].get("fileform").lower())
		if not overwrite_plots and os.path.exists(firstfilepath):
			verbPrint("Plot already exists. Skipping...")
			timeslider_status = TimeSliderNextState() # Advance to next state
			continue
		else:
			# file doesnt exist, so we should start the plotting here
			timeslider_status = TimeSliderNextState()
			break
		
		#check timeslider successfully executed
		if not timeslider_status:
			raise SystemError("TimeSlider could not be advance!")
		
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
			timeslider_status = TimeSliderNextState() # Advance to next state
			continue

		# if timeslider failed, return error
		if not timeslider_status:
			raise SystemError("TimeSlider could not advance!")

		print("Plotting file " + hdf5files[i])
		TimeSliderNextState()
		SaveWindowAtts = SaveWindowAttributes()
		SaveWindowAtts.outputToCurrentDirectory = 0
		SaveWindowAtts.outputDirectory = config["Output"]["output_plot_path"]
		SaveWindowAtts.fileName = str(variableToPlot) + ('%04d' % i)
		SaveWindowAtts.family = 0 # needs to be enforced again
		SetSaveWindowAttributes(SaveWindowAtts)	
		SaveWindow()
	
	# close the plots and check they were closed
	plotdelete_status = DeleteAllPlots()
	if not plotdelete_status:
		raise SystemError("Plot could not be deleted!")



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
