## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory 

import os, sys, shutil
import subprocess
import hashlib #for computing checksums

#add path to postprocessing module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..','..', 'PostProcessing')))
import Source

def download_dataset():
    """Check if the test datasets are downloaded. If not, download them from Zenodo.
    """

    #Zenodo ID for datasets
    #zenodo_id = 13257909
    zenodo_id = 13289221
    N_datasets = 12

    #download all test hdf5 data
    if not os.path.exists("./testdata"):
        os.mkdir("./testdata")
    for inx in range(N_datasets+1):
        padded_counter = str(inx).zfill(6)

        current_dataset_name = "test_dataset" + padded_counter + ".3d.hdf5"
        dataset_relative_path = "./testdata/" + current_dataset_name
        if os.path.isfile(dataset_relative_path):
            continue
        else:
            #download test hdf5 data. Only ~13 MB    
            zenodo_url = "https://zenodo.org/records/" + str(zenodo_id) + "/files/" + current_dataset_name
            wget_cmd = ["wget", "--quiet", "--no-check-certificate", "--output-document=" + dataset_relative_path, zenodo_url]
            zenodo_status = subprocess.run(wget_cmd)
            if not zenodo_status:
                print("Failed to download dataset!!")

    print("Finished downloading test hdf5 data!")

    if os.path.isfile("Integrals.dat"):
        print("Integral test data already downloaded!")
    else:
        #download test integrals data. only ~200 KB
        zenodo_url = "https://zenodo.org/records/" + str(zenodo_id) + "/files/Integrals.dat"
        wget_cmd = ["wget", "--quiet", "--no-check-certificate", "--output-document=" + "./testdata/Integrals.dat", zenodo_url]
        zenodo_status = subprocess.run(wget_cmd)
        if not zenodo_status:
            print("Failed to download test data! Exiting...")
            sys.exit("Failed to download test data!")

        print("Finished downloading test integral data!")
    
    ##Check data integrity with MD5 hash
    hdf5_hashes = ["1465595e89f6ae94a982a9ac92b9051f",
                    "bf4de829c8c843880702faef930b4ea3",
                    "f487ea79ef1c0964f0c766dd78b66d97",
                    "8b6880da47fce7e88bc86f866ad030af",
                    "f6a5d6828ef3134b3ff9991b122cb696",
                    "fad226297e44ebb40dd7cc1df94b7579",
                    "63108dc6e573e538e2eb66a7ee24f183",
                    "29f658e3d90e60bb873202e2a5f3a983",
                    "48a93962063fe7e4e288e0cdd77a97f5",
                    "22b8b0b11f56dbf02d17c29029321e6d",
                    "f203ba25284d7e6256a815a102506549",
                    "04dbc0e2684c508f48eebd524e62d78c",
                    "279872c5aa7ef9031a7a29d4a335af52"]
    dat_hash = "978a44cec209740c4f257afe893118c9"


    #generate checksums of downloaded data
    for inx in range(N_datasets+1):
        filename = "./testdata/test_dataset" + str(inx).zfill(6) + ".3d.hdf5"
        hdf5_checksum = hashlib.md5(open(filename, 'rb').read()).hexdigest() #compute md5 checksum

        #verify checksums against known values
        if not hdf5_checksum == hdf5_hashes[inx]:
            print("HDF5 test data integrity check failed!")
            sys.exit("HDF5 test data integrity check failed!")

    dat_checksum = hashlib.md5(open("./testdata/Integrals.dat", 'rb').read()).hexdigest()
    if not dat_checksum == dat_hash:
        print("Integral test data integrity check failed!")
        sys.exit("Integral test data integrity check failed!")
    
def LineTest():
    """Execute the test routines for the line plot postprocessing module

    Returns:
        int: success status of the test
    """

    #line number for the 'one_plot' parameter
    line_number = 8

    # First generate individual plots for each variable
    with open("IntegralsTestParams.ini", "r") as f:
        config_lines = f.readlines()
    
    config_lines[line_number] = "one_plot = 0\n"

    #write back to the .ini file
    with open("IntegralsTestParams.ini", "w") as f:
        f.writelines(config_lines)

    #now execute the test
    print("Executing single plot test...")
    subprocess_failure_status = subprocess.call(["python", "../../PostProcessing/IntegralsPlot.py", "IntegralsTestParams.ini"])
    if subprocess_failure_status:
        print("subprocess status code: ", subprocess_failure_status)
        raise RuntimeError("Line single plot test failed to execute!")

    #now execute the same test, but combine all plots into one by changing the 'one_plot' parameter
    config_lines[line_number] = "one_plot = 1\n"
    with open("IntegralsTestParams.ini", "w") as f:
        f.writelines(config_lines)

    print("Executing combined plot test...")
    subprocess_failure_status = subprocess.call(["python", "../../PostProcessing/IntegralsPlot.py", "IntegralsTestParams.ini"])
    if subprocess_failure_status:
        print("subprocess status code: ", subprocess_failure_status)
        raise RuntimeError("Line combined plot test failed to execute!")

    #compute the md5 checksum of the generated movie
    lineplot_checksum = hashlib.md5(open("./integral_plots/combined.png", 'rb').read()).hexdigest()

    #compute known checksum
    known_checksum = "3fd8d0836eac351c255516b46cf1feb2"
    

    #return both status codes
    return lineplot_checksum == known_checksum

def TwoDTest():
    """Executes the test routine for the 2D postprocessing module

    Returns:
        boolean: success status of the test
    """

    # Since this is a test script, we remove previous results 
    if os.path.exists("2d_plots"):
        print("Removing previous 2d_plots")
        shutil.rmtree("2d_plots")
    if os.path.exists("2d_movies"):
        print("Removing previous 2d_movies")   
        shutil.rmtree("2d_movies")  

    # Execute the tests
    print("Executing 2-D plot test...")
    subprocess_failure_status = subprocess.call(["visit", "-nowin", "-cli", "-s", "../../PostProcessing/VisitPlot.py", "2DPlotParams.ini"])
    if subprocess_failure_status:
        print("subprocess status code: ", subprocess_failure_status)
        raise RuntimeError("2-D plot test failed to execute!")

    #Compute the md5 checksum of the generated movie
    twodmovie_checksum = hashlib.md5(open("./2d_movies/Asquared.mp4", 'rb').read()).hexdigest() #compute md5 checksum
    known_checksum = "ee66ca9a2d8cb42d414c72a074d7fac1"
    twod_status = (twodmovie_checksum == known_checksum)

    return twod_status

def ThreeDTest():
    """Executes the test routine for the 3D postprocessing module

    Returns:
        boolean: success status of the test
    """
    print("\n\n\n\n")

    # Since this is a test script, we remove previous results
    if os.path.exists("3d_plots"):
        print("Removing previous 3d_plots")
        shutil.rmtree("3d_plots")
    if os.path.exists("3d_movies"):
        print("Removing previous 3d_movies")
        shutil.rmtree("3d_movies")

    # Execute the tests
    print("Executing 3-D plot test...")
    subprocess_failure_status = subprocess.call(["visit", "-nowin", "-cli", "-s", "../../PostProcessing/VisitPlot.py", "3DPlotParams.ini"])
    if subprocess_failure_status:
        print("subprocess status code: ", subprocess_failure_status)
        raise RuntimeError("3-D plot test failed to execute!")

    #Compute the md5 checksum of the generated movie
    threedmovie_checksum = hashlib.md5(open("./3d_movies/Asquared.mp4", 'rb').read()).hexdigest() #compute md5 checksum
    known_checksum = "3f9544bd1e076d7db1c3e623fffc8958"  
    threed_status = (threedmovie_checksum == known_checksum)                    

    return threed_status
    

if __name__ == "__main__": 

    #download test dataset from Zenodo
    download_dataset() 

    # Check if we're running VisIt
    try:
        __visit_script_file__
        running_visit = True
    except NameError:
        running_visit = False

    #if we're not running VisIt, then we're running python.
    if not running_visit:
        print("\n\n##############################")
        print("Running line test with python")
        print("##############################\n\n")
        line_success_status = LineTest()

        if line_success_status:
            print("\n\n###############")
            print("Tests passed!")
            print("###############\n\n")
        else:
            print("\n\n###############")
            print("Line test failed!")
            print("###############\n\n")

    elif running_visit:
        print("\n\n##############################")
        print("Running 2D and 3D test with VisIt")
        print("##############################\n\n")

        #run postprocessing for VisIt modules
        twod_success_status = TwoDTest()
        threed_success_status = ThreeDTest()

        success_status = (twod_success_status and threed_success_status)

        if success_status:
            print("\n\n###############")
            print("Tests passed!")
            print("###############\n\n")
        else:
            if not twod_success_status:
                print("\n\n###############")
                print("2D test failed!")
                print("###############\n\n")
            if not threed_success_status:
                print("\n\n###############")
                print("3D test failed!")
                print("###############\n\n")

        sys.exit(0 if success_status else 1)
    


