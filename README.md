# Bayesian_SSUR

### Dependencies
This software has been tested on UNIX systems only, on Windows it might work using Bash for Windows under the assumption that a recent version of the gcc compiler is available.
It depends heavily on Armadillo ( http://arma.sourceforge.net/ ) and on some boost libraryes ( https://boost.org ).
See the projects respective pages for installation.

## Usage
Run the sample_conformant_data.R script to simulate some synthetic data.
This will produce a bash script to call the program on those data as well; results will be presented in a 'results' subdirectory.

As prompted by the script, to run simply input

    make && chmod +x call.sh && ./call.sh
    
in the console. Run

    make clean
    
to cleanup. NOTE THIS WILL ERASE THE RESULTS AS WELL.

### Reference pre-print
https://www.biorxiv.org/content/early/2018/11/11/467019
