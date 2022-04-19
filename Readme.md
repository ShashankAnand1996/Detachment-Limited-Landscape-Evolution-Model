# Landscape Evolution Model (LEM)

- This Python code implements an algorithm to implictly and efficiently compute the erosion term for the detachment-limited landscape evolution model (LEM).
- The numerical results using this algorithm have been tested for accuracy based on analytical predictions in the publication referred below.
- The resulting numerical scheme allows us to achieve accurate steady-state solutions in conditions of high erosion rates leading to heavily dissected landscapes.

# Publication

The peer-reviewed research article can be found here: https://doi.org/10.1016/j.envsoft.2020.104804.

# Installation

- The code is presented as a *jupyter notebook*, which can be installed using using pip as `pip install notebook` or using Conda `conda install -c conda-forge notebook`.
- Required python libraries are *numpy* and *scipy* libraries, which can installed as `pip install numpy` and `pip install scipy`, respectively. 
- The presented algorithm is compatible with the *landlab* modeling environment in Python, which can be installed as `pip install landlab
` or using Conda `conda install landlab -c conda-forge`.
- The jupyter notebook also uses *display* module from the *ipython* library, which can be installed as `pip install ipython` or using Conda `conda install -c anaconda ipython`.

# Code structure

Each number written in the list below explains the operations and functions for the particular Cell of the Jupyter Notebook.
1. Import relevant packages from the above-mentioned libraries.
2. Import quadrature points and weight functions arrays (**s200.mat** and **w200.mat**), define Python functions for trial and test functions used in the weak formulation, and define Python functions to numerically obtain the slope and its derivatives for generic values of _m_ and _n_.
3. Select values of parameters (_C<sub>I</sub>,m,n,L,U,k_), global variables and arrays used in the solver.
4. For every Channelization Index (_C<sub>I</sub>_) and wavenumber (_k_), solve the eigen-value problem using spectral Galerkin technique with numerical quadrature.
5. Extract the critical Channelization Index _C<sub>I<sub>cr</sub></sub>_  and the corresponding fastest growing (positive) spatial frequency (_k<sub>cr</sub>_).
6. Save relevant arrays for plotting and further investigations.

# Results

* The "Results" folder in this github directory contains all the steady-state solutions shown in figures of the manuscript.
* Numpy arrays with keywords **acc** and **ele** indicate the total drainage area and elevation array respectively at the steady-state.
* Numpy arrays with keyword **diff** is the array of list with
  * diff[0] - time-step
  * diff[1] - maximum change in elevation at a node from previous time-step
  * diff[2] - change in mean elevation from previous time-step
