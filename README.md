sglake-parallel

[![DOI](https://zenodo.org/badge/428685309.svg)](https://zenodo.org/badge/latestdoi/428685309)

Authors: Aaron Stubblefield (Columbia University) with help/advising from Marc Spiegelman (Columbia University).

# Overview
This repository contains FEniCS python code for simulating subglacial lake shoreline
migration and free-surface evolution over time in two or three spatial dimensions. The model is
isothermal Stokes flow with nonlinear ("Glen's law") viscosity. The contact
conditions that determine whether ice remains in contact with the bed or
detaches are enforced with a penalty functional. The model runs in parallel with mpi4py (see below).

The default configuration is used to create synthetic data for the subglacial inversions in
the repository https://github.com/agstub/subglacial-inversion.

# Dependencies
## Required dependencies
As of this commit, this code runs with the latest FEniCS Docker image (https://fenicsproject.org/download/).
Docker may be obtained at: https://www.docker.com/. To run the Docker image:

`docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:current`


## Optional dependencies

1. FFmpeg (https://www.ffmpeg.org/) can be used, along with **process.py**,
to create a video of the evolving free surface geometry over time. 

2. ParaView is useful for visualizing velocity/pressure solutions to the Stokes equations (https://www.paraview.org/).

# Contents

## 1. Source files
The model is organized in a series of python files in the *source* directory as follows.

1. **geometry.py** contains the geometric description of the bed and initial ice-water interface.

2. **params.py** contains all of the model parameters and model options.

3. **stokes.py** contains the Stokes system solver and related functions.

4. **mesh_fcns.py** contains functions that solve the surface kinematic equations and move the mesh.

5. **boundary_conds.py** contains functions that mark the mesh boundary and apply boundary conditions.

6. **hydrology.py** contains the subglacial lake volume change
timeseries.

7. **main.py** runs the model. It contains the time-stepping loop that
calls the Stokes solver and mesh-related functions at each timestep, and saves the output.

8. **post_process.py** and **realtime_process.py** extract elevation and basal vertical velocity fields
from the FEniCS velocity solution and mesh. The only difference is that **realtime_process.py**
does this while the simulation is running, rather than afterwards.

9. **plotting.py** creates png images of the elevation anomaly, basal water thickness,
and basal vertical velocity at each time step. This is called by
**post_process.py**.



# Running the code
To run the code:

1. Run the FEniCS Docker image.

2. In Docker, run the main file from the parent directory: `python3 ./source/main.py`

The code can be run in parallel
via `mpirun -np num_proc python3 ./source/main.py` where `num_proc` is the
number of processes.

# Output

Model output is saved in a *results* directory. This includes

1. the Stokes solution vtk files (*results/vtkfiles* subdirectory) (if turned on in params.py),

2. hdf5 files (*results/h5files* subdirectory) that are used in **post_process.py** to
plot quantities of interest.

# Creating synthetic data for inversions
After the code has run to completion in the default configuration,
just run `python3 ./scripts/process.py` to create numpy arrays
`h_true.npy` and `wb_true.npy`, along with a figure of the solution at
a few time steps (fig8.png) and a collection of png's of the solution
at every time step (saved in a directory called *movie*).
These numpy arrays are then used as example synthetic data in the subglacial inversion repository
https://github.com/agstub/subglacial-inversion.

# Model options
Model options and parameters are set in the **params.py** file.
