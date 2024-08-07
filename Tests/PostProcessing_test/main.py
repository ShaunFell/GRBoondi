## GRBoondi 2024
## Copyright 2024, Shaun Fell
## Please refer to LICENSE in GRBoondi's root directory

import os, sys
import subprocess


def main():
    """This script tests the PostProcessing module.
    """
    success_status = 1

    #add path to postprocessing module
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..','..', 'PostProcessing')))
    
    #download test hdf5 data. Only ~13 MB    
    zenodo_url = "https://zenodo.org/records/13247864/files/GeneralizedProcap_000000.3d.hdf5"
    zenodo_status = subprocess.run(["wget", "--quiet", "--no-check-certificate", "--output-document=test_dataset.3d.hdf5", zenodo_url]) 
    if not zenodo_status:
        print("Failed to download test data! Exiting...")
        sys.exit("Failed to download test data!")

    #run postprocessing




    return success_status

if __name__ == "__main__":
	main()