# sfec13d makefile

#compiler = /opt/local/bin/gfortran
compiler = /usr/local/bin/gfortran-12
options = 
output = -o sfec
programs: nrtype.f90 nrutil.f90 geometry.f90 modules.f90 sfec13d.f90 
	$(compiler) $(options) $(output) nrtype.f90 nrutil.f90 geometry.f90 modules.f90 sfec13d.f90

clean: 
	rm -f sfec *.mod


