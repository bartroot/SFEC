# SFEC

The modified Spectral Finite Element Code (SFEC) can be linked with the GSH repository, using this modified version of SFEC.

## Requirements

The code need to be compiled with a gfortran compiler. The newest version is advised. You can specify your compiler in the makefile.

## Structure and usage

In this repository you can find:

Initial Fortran codes:
- `sfec13d.f90`: Main source code that contains the Spectral Finite Element Stokes and Poisson solver
- `geometry.f90`: Specifies the geometry of the model
- `modules.f90`: Fortran numerical modules that are used in the main code
- `nrtype.f90`: Type definition file used in the main code
- `nutil.f90`: Utility file used in the main code
- `makefile`: makefile to compile the fortran code
- `sfec.inp`: input file used by SFEC to run through different scenarios

Scripts that are needed to link to GSH code using the 'Model' format to describe the density model:

- `construct_SFEC_inp.m`: Initialise the SFEC.inp file that is the input file for the fortran SFEC code
- `construct_geometry_mod.m`: rewrite the geometry fortran file, such that SFEC can be used for any planetary density shell
- `create_tom_files.m`: create density tomography files that are the driving force of the convection code. This code creates `outputmax.tom` file that is used in sfec.inp to describe the density forcing in the Stokes solution.

To compile the code, go to the directory and type

`./make`

Several .mod files are generated and an executable:

`./sfec`: This runs the whole code.

When `./sfec` is run, it generates two output files, containing the gravity and dynamic topography signatures of the convecting shell:

- `geoid.dat`
- `geoid_coeff.dat`

### Modelling scenarios

Mars: Mantle plume dynamics linked to flexure theory derived lithosphere model 

### Keep in Mind

For now this package can be run stand alone or linked with the GSH code for full planetary gravity field modeling.

## Authors (Standing on the shoulders of giants!)

This Software has been developed on ideas and software from the following developers/contributors:

- **Nicola Tosi** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)[0000-0002-4912-2848](https://orcid.org/0000-0002-4912-2848), DLR Berlin (Developer)
- **Zdenek Martinec** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)[0000-0002-7036-2571] (https://orcid.org/0000-0002-7036-2571), Dublin Institute for Advanced Science (Developer)
- **Bart Root** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)[0000-0001-7742-1434](https://orcid.org/0000-0001-7742-1434), Technische Universiteit Delft (Contributor)
- **Youandi van der Tang**, Technische Universiteit Delft (Contributor)
- **Maxime Vincent**, University of Paris (Contributor)  

## License

The contents of this repository are licensed under a **GNU General Public License v3.0** (see [LICENSE](https://github.com/bartroot/GSH/blob/main/LICENSE.md) file).

Copyright notice:

Technische Universiteit Delft hereby disclaims all copyright interest in the program “SFEC”. SFEC is a Fortran/MATLAB package to do perform the Stokes solve of a convective density shell, written by the Author(s). 
Henri Werij, Faculty of Aerospace Engineering, Technische Universiteit Delft. 

© 2021, B.C. Root

## References

If you would like to reuse the code, please cite the following:

Tosi, N. (2007), Numerical modeling of present-day mantle convection, PhD thesis, Charles University, Prague, (https://geo.mff.cuni.cz/theses/2008-Tosi-PhD.pdf)

or

Tosi N., Z. Martinec (2007). Semi-analytical solution for viscous Stokes flow in two eccentrically nested spheres. Geophysical Journal International, 170, 1015-1030, doi:10.111/j.1365-246X.2007.03482.x

## Would you like to contribute?

If you have any comments, feedback, or recommendations, feel free to **reach out** by sending an email to b.c.root@tudelft.nl

If you would like to contribute directly, you are welcome to **fork** this repository.

Thank you and enjoy!
