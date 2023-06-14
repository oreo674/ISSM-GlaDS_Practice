
LIBS=SaveNetcdf \
     USF_CouplingSSA

NETCDF=$(NETCDFHOME)
F90=gfortran
PROG=../PROG

all: $(LIBS)
    

SaveNetcdf: $(PROG)/SaveNetcdf_SHMIP.F90
	elmerf90 -I$(NETCDF)/include  $^ -o $@ -L$(NETCDF)/lib -lnetcdff -lnetcdf
       
USF_CouplingSSA: $(GlaDS)/USF_CouplingSSA.f90 
	elmerf90  $^ -o $@ 

Mesh:
	gmsh -1 -2 mesh_unst.geo
	ElmerGrid 14 2 mesh_unst.msh -autoclean

clean:
	rm -f $(LIBS)
