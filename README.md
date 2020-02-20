# Scripts for computational chemistry

Here is presented a couple of scripts intended to facilitate the data treatment of computational chemsitry calculations. Some programs support results for different softwares, others are package specific. Each script is intended to be self contained, this may cause bits of code being repeated in different programs, but removes dependence among them. For more details, check each program's help.

The available scripts are:
- [scan](./doc/scan.md): Takes the vertical excitation energies of different geometries following a specific parameter (bond distance, angle, dihedral) or no parameter at all (like for a breathing mode). Currently support only for ORCA.
- [spectrum](./doc/spectrum.md): Reads the vertical excitation energies from a TD-DFT calculation and calculates gaussian broadening to reproduce experimental result. Currently support for ORCA and Gaussian09.
- [orca extractor](./doc/orca_extractor.md): Reads ORCA output files and returns the asked results in a friendly output.

The scripts are being written as the necessity arrives and the time allows.

## Usage

There is no need for instalation, just copy the script in any desired folder and call it to use. For a more detailed explanetion on using the scripts, look the script specific exlpanation in the [doc](./doc) files.

## Credits

I'm a chemsitry student at UFRJ pursuing a Master's degree. This project is intended to help me practice python, managing a repository and facilitate the lab work. Any tips and suggestions on what to add to the scripts are welcome.