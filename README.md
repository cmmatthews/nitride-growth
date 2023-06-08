# nitride-growth
A Kinetic Model for Binary &amp; Ternary III-Nitride Epitaxy
* This code implements a kinetic model for the growth of III-nitride semiconductor thin films. The current version assumes that the mechanisms are independent of temperature, making the current form of this model most appropriate for the growth of AlGaN. Additionally, the model can simulate the growth of binary III-nitrides (again assuming no thermal effects).
* There are two main use cases for this model: 1) investigate the physical mechanisms that control III-nitride epitaxy, and 2) simulate the growth of III-nitride thin films/structures.
* 30 seconds of simulated growth typically runs in 1.5-2 seconds on a machine with a 3.6 GHz processor and 16 GB of RAM 

## Quickstart
1. Run ``modelAthermal.m`` in MATLAB
   1. For any run with at least one of the variables (``fluxN``, ``IIIV``, ``xA``, ``rPMLC``, or ``rLCDr``) defined as a vector, the script will generate and solve multiple growth scenarios
   1. Runs with only scalar variables will only generate and solve one growth scenario
1. The output figures will be saved as both .fig and .tif files in ``[Model Path]/output/[DateTime]/``
1. Raw data for the run will be in the ``dataTables`` variable in the MATLAB workspace

## Citing

If you use this model in an academic context, please cite following publication(s):

``
@article{,
  author    = {},
  title     = {},
  journal   = {},
  volume    = {},
  year      = {},
  url       = {}
}
``

## Running Custom Simulations

### Changing Physical Mechanisms
1. Open ``modelAthermal.m`` in MATLAB
1. Change rate coefficients on lines 59-63 (in Physics Variables Section)
   1. Note that each VCS coefficient (``rPMLC`` and ``rLCDr``) can be a scalar or vector to enable solving multiple growth scenarios at once
1. Run ``modelAthermal.m`` in MATLAB
1. Output will be saved as described in the Quickstart section

### Changing Growth Varialbes
1. Open ``modelAthermal.m`` in MATLAB
1. Change growth variables on lines 12-44 (in User Variables Section)
   1. Note that ``fluxN``, ``IIIV``, and ``xA`` can all be scalar or vectors to enable solving multiple growth scenarios at once
   1. The variable ``tSubs`` can also be scalar or a vector, but is currently unused in this version of the model
1. Run ``modelAthermal.m`` in MATLAB
1. Output will be saved as described in the Quickstart section

### Changing Output Settings (Not Recommended)
1. Open ``modelAthermal.m`` in MATLAB
1. Change output variables on lines 48-54 (in User Variables Section)
   1. **All of these output options are legacy code and disabled by default, thus enabling them may have unexpected results**
   1. ``outputSettings`` is used to enable the output of the Runge-Kutta solution to the console, a file, both or neither
   1. ``fileOut`` has a similar effect to setting ``outputSettings = 'file'``, but is implemented in a different place in the code
   1. ``genPPT`` generates a PowerPoint slide with some of the figures that are created as part of the normal output. This uses the slide templates saved in ``[Model Path]/PowerPoint Templates``
1. Run ``modelAthermal.m`` in MATLAB
1. Output will also still be saved as described in the Quickstart section

## Basic Program Flow
1. Defining the growth scenarios (stored in ``system``)
    1. An empty 5-D (or 6-D) matrix is created based on the values of ``fluxN``, ``IIIV``, ``xA``, ``rPMLC``, and ``rLCDr`` (plus ``tSubs`` although it is unused) to store the results of simulating growths in each combination of these variables. (line 134 of ``modelAthermal.m`` -> ``genSystem()``)
    1. Within the `genSystem` function, the rate equations described in the paper associated with this model are implemented as a cell array of anonymous functions which each call the relevant rate functions in the ``Rate Functions`` section of ``modelAthermal.m`` (lines 736-925)
    1. Each of the rate functions also calls relevant functions for each physical rate modeled here, which are implemented in the ``Rate Definitions`` section of ``modelAthermal.m`` (starting on line 925)
1. Solving each growth scenario
    1. All growth scenarios are evaluated one at a time in the loop structure starting on line 168 of ``modelAthermal.m``
    1. Each scenario's system of rate equations is solved by the Runge-Kutta solver (``RK4()``), and the output is stored in the 5-D (or 6-D) ``dataTables`` variable.
        1. The ``RK4`` function also supports a callback method that evaluates after every timestep in the main RK4 loop, which is currently used to update the semiconductor crystal composition profile as the system is being solved. *This is likely not important for the current iteration of this model, but may be crucial when thermal rates are incorporated*
1. Plotting the results
    1. All of the simulated growths are plotted one at a time in the loop structure starting on line 252 of ``modelAthermal.m``
    1. The ``plotMLs`` function is used to generate a 2x2 plot of the # of atoms in each layer of the structure (semiconductor film/crystal, PM-ML, LC-ML, and droplets adlayer) vs. time, plus a second plot isolating the PM-ML evolution
    1. The ``plotCrystal`` function is used to generate two histograms that show the composition of each monolayer of semiconductor crystal that was grown during the simulated growth, one with labels and one without.
    1. The four plots generated in this section are saved as .fig and .tif to the ``output`` directory by default.
1. End of script error checking
    1. Lines 338-348 of ``modelAthermal.m`` are used to compare the input and output values for total # of atoms in the system, to check that mass-balance is met.

## III-Nitride Growth Model File Descriptions

### Base Files
``modelAthermal.m``
* main script
* models III-nitride growth for a set of user defined growth variables
* outputs plots for: # of atoms vs. t (grown crystal and 3 adlayers) and layer-by-layer composition profile
* raw data can also be saved manually or using various output flag variables

``myStep.m``
* function definition for making step functions for flux profiles

``genMatParams.m``
* function definition to define lattice parameters for III-nitride binaries
* saves output to materialParams.mat

``materialParams.mat``
* data file storing III-nitride binary lattice parameters
* created/updated by genMatParams.m

### Utility Files

``utils/barExtract.m``
* script to be used to extract layer-by-layer data from a set of save .fig files
* this data is used for layerStats.m

``utils/layerStats.m``
* script to calculate average composition and standard dev of composition in super lattice sublayers, which are defined in the first few lines of the script
* this data is used for optimizationPlots.m

``utils/optimizationPlots.m``
* script to plot super lattice layer statistics as a function of the model's free parameters
* colors of data points represent values of each statistic

### PowerPoint Templates

``PowerPoint Templates/...``
* contains blank PPT files for one form of output that's built in to the main script
* untested for recent versions, may be deprecated

### Output Files

``output/...``
* stores output .fig and .tif files of main script
* created by modelAthermal.m if it doesn't exist

## License
-----

This model is released under the BSD license, reproduced in the LICENSE file in this directory.
