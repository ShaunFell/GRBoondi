## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

import os, glob, sys
import Source
import configparser

## Note:
# we must put user-installed modules inside the functions that use them since VisIt doesn't support additional modules
# and it's a pain to set it up for that. For simplicity, we put them there.

## flag for VisIt
try:
    import visit
    __visit_imported = True
except:
    __visit_imported = False

## flag for Python
try:
    import matplotlib
    import pandas
    __extra_python = True
except:
    __extra_python = False

#decorator for functions that require VisIt
def require_visit(fn):
    def wrapper(*args, **kwargs):
        if not __visit_imported:
            raise RuntimeError("VisIt not found. Aborting.")
        import visit
        return fn(*args, **kwargs)
    return wrapper

#decorator for functions that require non-standard python packages, e.g. matplotlib and pandas
def require_python(fn):
    def wrapper(*args, **kwargs):
        if not __extra_python:
            raise RuntimeError("Non-standard python packages not found. Aborting.")
        import matplotlib.pyplot
        import pandas
        return fn(*args, **kwargs)
    return wrapper


def get_hdf5_file_list(config):
    """Find all the hdf5 files and return them as list of absolute path strings

    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters

    Raises:
        FileNotFoundError: Error raised if hdf5 path isn't found
        FileNotFoundError: Error raised if no plot files were found

    Returns:
        list: list of plot file names
    """
    hdf5_path = config["Header"]["hdf5_path"]
    filenames = config["Header"]["plot_header"] + "*.hdf5"

    if not os.path.exists(hdf5_path):
        raise FileNotFoundError("HDF5 path not found: " + hdf5_path)

    #get list of file using regex 
    hdf5files = glob.glob(os.path.join(hdf5_path, filenames))
    hdf5files.sort() # sort the files by number

    #check if files exist
    if len(hdf5files) == 0:
        raise FileNotFoundError("No HDF5 files found in " + config["Header"]["hdf5_path"])

    return hdf5files

def get_config(arg_list):
    """safely get config object from command line arguments

    Args:
        arg_list (list): list of arguments from command line

    Raises:
        RuntimeError: please specify a parameter file

    Returns:
        configparser.ConfigParser: instance of a ConfigParser class that holds the users parameters
    """

    if len(arg_list) < 2:
        raise RuntimeError("Please specify a parameter file")
    
    config = configparser.ConfigParser()
    config.read(arg_list[1])

    if __visit_imported: #If we're running VisIt, then these are the parameter sections we expect
        mandatory_visit_headers = ["Header", "EngineConfig", "VariableData", 
                            "AnnotationConfig", "PlotConfig", 
                            "ViewConfig", "Output"]
    elif __extra_python: #If we're running via python, then these are the parameter sections expect
        mandatory_visit_headers = ["Header", "VariableData", "PlotConfig", "Output"]


    loaded_config_sections = config.sections()  

    for header in mandatory_visit_headers:
        if header not in loaded_config_sections:
            raise configparser.NoSectionError("Parameter file error. Missing section: " + header)

    return config

def create_output_dirs(config):
    """create output directories if they don't exist

    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters
    """
    plot_path = config["Output"].get("output_plot_path", fallback = "./plots")
    movie_path = config["Output"].get("output_movie_path", fallback = "./movies")

    if not os.path.exists(plot_path):
        os.mkdir(plot_path)
    if not os.path.exists(movie_path):
        os.mkdir(movie_path)

def plot_func_selector(type):
    """selector function for type of plot to generate

    Args:
        type (str): string defining type of plot function to choose. Options: '2d' and '3d'

    Raises:
        ValueError: Invalid plot type

    Returns:
        func: appropriate plotting function
    """

    if type == '2d':
        return Source.TwoD.generate_2dslice_plot
    elif type == '3d':
        return Source.ThreeD.generate_3dvolume_plot
    else:
        raise ValueError("Invalid plot type. Must be either '2d' or '3d'")

def PlotFiles(config):
    """	Find all the plot files and return them as list of absolute path strings

    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters

    Raises:
        FileNotFoundError: Error raised if no plot files were found

    Returns:
        list: list of plot files
    """

    #path to hdf5 files plus plot file header string
    filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])

    #get list of file using regex
    files = glob.glob(filename_prefix+ "*.3d.hdf5")

    #filter the file list to ensure regex-captured files are actual plot files
    plot_files = [x for x in files if config["Header"]["plot_header"] in x]

    #ensure plot files exist in provided directory
    if len(plot_files) == 0:
        raise FileNotFoundError("No Plot Files Found!")

    return plot_files

def MultipleDatabase(config):
    """Determine whether multiple plot files need to be opened

    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters

    Returns:
        boolean: boolean flag determining if multiple plot files need to be opened
    """

    if len(PlotFiles(config))>1:
        print("Opening multiple plot files")
        return True
    else:
        print("Opening single plot file")
        return False

@require_visit
def Open_Database(config):
    """Opens a visit database object

    Args:
        config (configparser.ConfigParser): instance of a ConfigParser class that holds the users parameters

    Raises:
        IOError: Could not open database!

    Returns:
        int: status code of opening the database
    """

    if MultipleDatabase(config):
        filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
        database_status = visit.OpenDatabase(filename_prefix + "*" + ".3d.hdf5 database", 0)
    else:
        database_status = visit.OpenDatabase(PlotFiles(config)[0], 0)

    if not database_status:
        raise IOError("Could not open database!")

    return database_status

@require_visit
def Close_Database(config):
    
    if MultipleDatabase(config):
        filename_prefix = os.path.join(config["Header"]["hdf5_path"], config["Header"]["plot_header"])
        database_status = visit.CloseDatabase(filename_prefix + "*" + ".3d.hdf5 database")
    else:
        database_status = visit.CloseDatabase(PlotFiles(config)[0])

    if not database_status:
        raise IOError("Could not close database!")

    return database_status

@require_python
def get_data_dataframe(filepath):
    """
    Read in the data file and return a pandas dataframe

    Args:
        filepath (str): absolute path to the data file

    Returns:
        pd.DataFrame: a pandas dataframe
    """

    #get the header lines
    #from GRBoondi, they should always start with a '#'
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("#"):
                header = line
            else:
                break #stop the loop if theres no more header liens

    header = header[1:].strip().split()
    
    #import the data to a pandas dataframe and return it. Includes error checking
    return pandas.read_csv(filepath, sep='\s+', names=header, engine="python", header=0, index_col=False)

@require_python
def setup_pyplot(config):
    """
    Setup the plot style
    """

    # turn off interactive plotting and set plot style
    matplotlib.pyplot.ioff()
    matplotlib.pyplot.style.use(config["PlotConfig"].get("plotstyle", "default"))

@require_python
def setup_pyplot_figure(): 
    """
    Setup the figure environment
    """

    matplotlib.pyplot.figure()
    matplotlib.pyplot.tight_layout()

@require_python
def select_plot_func(scaling):
    """
    Select the plot function based on the scaling

    Args:
        scaling (str): "linear", "loglog", "loglinear", "linearlog"

    Returns:
        function: the plot function from pyplot
    """

    if scaling == "loglog":
        return matplotlib.pyplot.loglog
    elif scaling == "loglinear":
        return matplotlib.pyplot.semilogx
    elif scaling == "linearlog":
        return matplotlib.pyplot.semilogy
    else:
        return matplotlib.pyplot.plot

@require_python
def save_pyplot_fig(config, varname):
    """
    Save the pyplot figure

    Args:
        path (str): absolute path to the output directory
        varname (str): name of the variable
    """

    path = config["Output"].get("output_plot_path", fallback = "./plots")
    save_format = config["Output"].get("fileform", fallback = "png").lower()
    save_dpi = config["Output"].getint("dpi", fallback = 600)
    save_filename = varname+".{}".format(save_format)

    #save to disk
    matplotlib.pyplot.savefig(os.path.join(path,  save_filename), format = save_format, dpi = save_dpi)


def make_movie(config, plotvariable):

    make_movie = config["Output"].getint("make_movie", fallback = 0)

    if MultipleDatabase(config) & make_movie:
        print ("Making a movie...")
        framerate = config["Output"].getint("movie_framerate", fallback = 5)
        verbosity = config["Output"].get("verbosity", fallback = "0")
        absolute_input_files = os.path.join(config["Output"]["output_plot_path"], plotvariable+"%04d.mp4")
        options = ' -vcodec libx264 -crf 25 -pix_fmt yuv420p '
        absolute_output_file = os.path.join(config["Output"]["output_movie_path"], plotvariable+".mp4")

        ffmpeg_cmd = 'ffmpeg -r ' + str(framerate) + ' -y -v ' + str(verbosity) + ' -s 1920x1080 -i ' + absolute_input_files + options + absolute_output_file
        ffmpeg_system_status = os.system(ffmpeg_cmd)
        
        #Check ffmpeg status to ensure command executed successfully
        if not ffmpeg_system_status == 0:
            raise 	OSError("ffmpeg failed. Could not make the movie. cmd: {0}".format(ffmpeg_cmd))


class VerbosityPrint:
    """
    class that allows printing using user-specified verbosity levels
    
    Instantiate class using a given verbosity level, then call the instance 
    by passing objects to be printed as arguments.

    Example:

        verbPrint = VerbosityPrint(1)
        verbPrint("I'm an example!")

        output: "I'm an example!"

        verbPrint = VerbosityPrint(0)
        verbPrint("I'm another example!")

        output: ""
    """

    def __init__(self, verbosity):
        self.verbosity = verbosity

    def __call__(self, *objects):
        if self.verbosity:
            print(*objects)
