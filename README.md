# optimal_sensing_FEAwing

This repository accompanies the paper "Nonuniform structural properties of wings confer sensing advantages" by Weber*, Babaei*, Mamo, Brunton, Daniel, & Bergbreiter (*equal contributions)

To run this code:
1) Download example datasets (or create your own).
2) Preprocess data files output from COMSOL.
3) Run sensor optimization code.

## 1) Download example data

Find example datasets (.csv files) at:
[https://doi.org/10.5061/dryad.fxpnvx0wq](https://doi.org/10.5061/dryad.fxpnvx0wq)

At the above link, you will also find a file which can be used to create datasets in COMSOL Multiphysics.  The data accompanying the paper was created using version 5.6.


## 2) Preprocess data files

Run the script `preprocess_comsol_data` to properly sort and trim the example datasets. 
Note that the code assumes the working directory is this Github repository (optimal_sensing_FEAwing) and that the example datasets are located directly within this directory.  You may need to change paths to data files in the code, particularly if you are using your own dataset.


## 3) Run sensor optimization code

Run the script `wing_sensors_comsol_main` to determine optimal sensor placement and classification accuracy for the desired dataset(s). Again, you may need to change paths for the code to properly run.

Note that this code relies on the `cvx` package, which can be found at [http://cvxr.com/cvx/](http://cvxr.com/cvx/).