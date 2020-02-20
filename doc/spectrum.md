# Spectrum

Script to read from the result files with the vertical excitation energies, and calculate the broadening for comparison of theoretical values with experimental UV/Vis data.

This script began with modifications of the [g09_spectrum.py](https://github.com/mdommett/compchem-scripts/blob/master/g09_spectrum.py) to work with ORCA as well, and thus began this project. The broadening gaussianis calculated following the [Gaussian instructions](http://gaussian.com/uvvisplot/).

## Usage
The only required arguments are the calculation result files, which in the moment can be ORCA's "\*.out" or Gaussian09's "\*.log". The script identifies automatically whether it is a ORCA or Gaussian output file. Different result files can be given simultaneously, returning a plot with one line per file.<br>
(spectrum&#46;py OUTPUT1 OUTPUT2 ...)

The program can also be specified manually with the respective keyword.
- prog: the supported options are "orca" and "gaussian".<br>
(spectrum&#46;py OUTPUT -prog orca|gaussian)

#### Controling the output
There are many optional arguments to be used. If none is used, the program will return nothing. All options can be used separatelly or in conjunction.

- sd: sets the standard deviation in eV used to calculate the broadening. The default value is 0.4 eV;<br>
(spectrum&#46;py OUTPUT -sd 0.01)
- rng: specifies the min and max value for the spectra (in nm);<br>
(spectrum&#46;py OUTPUT -rng 200 900)
- sticks: used in the gnu and matplot plots. Enables plotting the spectrum with vertical bars corresponding to the theorefical excitations, which were used to calculate the spectrum.<br>
(spectrum&#46;py OUTPUT -sticks)

#### Plotting
These options are intended for fast visualization, they'll show the plot but won't save it.

- gnu: creates a plot using gnuplot.<br>
(spectrum&#46;py OUTPUT -gnu)
- mpl: creates a plot using matplotlib.<br>
(spectrum&#46;py OUTPUT -mpl)

#### Saving the result
These options are used to save the spectrum in the preferred format. All options can be used with or without argument. If the option is called with no argument, the default name used for the spectrum is "spectrum" (I know, not that creative).

- save: saves the matplot plot in a pdf format. The "-mpl" option don't have to be used to be able to generate and save the pdf;<br>
(spectrum&#46;py OUTPUT -save)
- dat: saves the data as a text file, separeted in columns. It can be read by another program to plot according to the user's need;<br>
(spectrum&#46;py OUTPUT -dat)
- csv: similar to saving with dat, but the resulting file is formated as a "comma separated values".<br>
(spectrum&#46;py OUTPUT -csv)

## TODO
#### Broadening
- Option to use Lorentzian broadening (is convolution hard to do? Let's find out...);
- Add options for the resulting spectrum unit.

#### Supported programs
- TURBOMOLE;
- Implementation of the [cclib](https://cclib.github.io/) package, enabling work with all cclib supported programs;

#### Outputting
- Hability to save in different formats.

#### Fixes
- Plotting spectra of multiple files with gnuplot and data files.