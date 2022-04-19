# Landscape Evolution Model (LEM)
- This Python code implements an algorithm to implictly and efficiently compute the erosion term for the detachment-limited landscape evolution model (LEM).
- The numerical results using this algorithm have been tested for accuracy based on analytical predictions in the publication referred below.
- The resulting numerical scheme allows us to achieve accurate steady-state solutions in conditions of high erosion rates leading to heavily dissected landscapes.

# Publication
The peer-reviewed research article can be found here: https://doi.org/10.1016/j.envsoft.2020.104804.

# Installation
- The code is presented as a jupyter notebook, which can be installed using using pip as `pip install notebook` or using Conda `conda install -c conda-forge notebook`.
- Required python libraries are numpy and scipy libraries. 
- The presented algorithm is compatible with the Landlab modeling environment in Python, which can be installed as .
- The jupyter notebook 


#Results

- The "Results" folder in this github directory contains all the steady-state solutions shown in figures of the manuscript.
- Numpy arrays with keyword **acc** and **ele** indicate the accumulation array and elevation array respectively at steady state.
- Numpy arrays with keyword **diff** is the array of list with
-- diff[0] - time-step
-- diff[1] - maximum change in elevation at a node from previous time-step
-- diff[2] - change in mean elevation from previous time-step
