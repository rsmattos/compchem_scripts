# Scan

Given a geometrical parameter and a set of output files, the script reads the energies in earch of them, the value of the parameter being evaluated and generate a plot or datafile for the energies. For now the only program supported is ORCA, tested with version 4.2.1.

## Usage

#### Parameter
The only required argument is the parameter being evaluated. The atoms involved have to be declared by its index of the cartesian coordinates starting with the first atom being counted as 1( not 0). The script only accepts one parameter argument at a time.

For the time being, the supported parameters are:

- bond distance: given any two atoms, calculates the distances along the scan. The atoms don't have to be actually bonded.<br>
(scan&#46;py -b|--bond ATOM1 ATOM2);
- angle: given any three atoms, calculates the angle. Again, the atoms don't have to be bonded. The order they are passed is importante, the angle calculated is in the ATOM2.<br>
(scan&#46;py -a|--angle ATOM1 ATOM2 ATOM3)
- dihedral: given any four atoms, calculates the dihedral angle they form. Once again, bonding doesn't matter, but ordering does. Calculates the angle between the planes characerized by ATOMS 1 2 3 and ATOMS 2 3 4.<br>
(scan&#46;py -d|--dihedral ATOM1 ATOM2 ATOM3 ATOM4)

#### Reading the outputs
There are different ways to tell the script what are the output files. By default, if nothing is specified, the script will look in the folder it's being run and it's subfolders for files with the extension "\*.out". This behaviour can be modified.

The script can be given a file, a list of files, a directory and a list of directories, as positional arguments. It can even be given some files and some directories. It will read the files and will look in the directories for output files with the given extension. The paths can be given as absolute or relative paths.

- To search in the current directory;<br>
(scan&#46;py -b 2 4);
- To search in a given directory;<br>
(scan&#46;py /path/to/directory -b 2 4).
Notice that the position doesn matter;<br>
(scan&#46;py -b 2 4 /path/to/directory)
- To read a given file, or a list of files.<br>
(scan&#46;py -b 2 4 /file/1 /file/2 /file/3 ...)

If the user preffers, it's also possible to create a file with the paths to the outputs, a file named OUTPUTS which contains the lines:

/path/to/output/1<br>
/path/to/output/2<br>
/path/to/output/3

This can be passed to the script as (scan -b 2 4 -p OUTPUTS). "-p" stands for short of "--paths_file".

The user is also given the option to change the extention. Using (scan&#46;py -b 2 4 -e .inp) will look for the files ending in "\*.inp" instead.

## TODO
#### Reading
- Add option to not look in the subfolders
#### Parameter
- Support for coordinates systems other than cartesian;
- Follow atom moving out of plane;
- General parameter, useful for vibrational modes, where the whole molecule is varying;

#### Supported programs
- Gaussian 09, since the data extraction is quite close to ORCA;
- TURBOMOLE,
- Implementation of the [cclib](https://cclib.github.io/) package, enabling work with all cclib supported programs;

#### Outputting
- Improve plotting with matplotlib;
- Create option to plot with gnuplot;
- Create option to output save a data file and csv;
- Create option to save plot as pdf and other formats.