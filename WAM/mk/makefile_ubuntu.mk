
# ubuntu on local pc:
FC = /disk1/local/openmpi_gfortran/bin/mpif90
FFLAGS = -O3 -ffree-form -fno-range-check -fno-second-underscore
INC = -I/disk1/local/openmpi_gfortran/include
NCDIR=/opt/netcdf-fortran-4.2

# vilje:
# FC = ifort
#flags='-free -c -O3 -xAVX -assume cc_omp  -traceback -ip -align all -ftz -fno-alias -no-prec-div -no-prec-sqrt -assume  dummy_aliases' 


#vpath %.f90 ../src/chief
#.SUFFIXES:
#.SUFFIXES: .f90 .o .mod


# rule for creating object and module files:
%.o %.mod: %.f90
	$(FC) $(INC) $(FFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(INC) $(FFLAGS) -c -o $@ $<

all: preproc wam pnetcdf

preproc: wam_special_module.o wam_model_module.o wam_mpi_module.o wam_file_module.o \
	wam_general_module.o wam_timopt_module.o wam_fre_dir_module.o \
	wam_interface_module.o wam_grid_module.o wam_boundary_module.o preproc_module.o \
	preproc_user_module.o wam_nest_module.o read_boundary_input.o wam_output_set_up_module.o \
	preproc.o read_topography.o read_preproc_user.o wam_mpi_comp_module.o wam_coordinate_module.o
	$(FC) $(INC) $(FFLAGS) $^ -L/usr/lib/openmpi/lib -lmpi -o $@

wam: wam_file_module.o wam_general_module.o wam_timopt_module.o wam_fre_dir_module.o \
	wam_interface_module.o wam_grid_module.o wam_current_module.o wam_model_module.o \
	wam_ice_module.o wam_output_module.o wam_wind_module.o wam_boundary_module.o \
	wam_source_module.o wam_propagation_module.o preproc_module.o wam_coldstart_module.o \
	wam_restart_module.o wam_initial_module.o wam_mpi_module.o wam_output_set_up_module.o \
	wam_topo_module.o wam_radiation_module.o wam_nest_module.o wam_user_module.o \
	wam_special_module.o read_topo_input.o wavemdl.o initmdl.o read_wam_user.o \
	print_wam_status.o read_wind_input_metno.o read_current_input.o wamodel.o read_boundary_input.o \
	read_ice_input_metno.o jafu.o wam_mpi_comp_module.o wam_assi_set_up_module.o \
	wam_assi_module.o wam_coordinate_module.o readsat.o chief.o
	$(FC) $(INC) $(FFLAGS) $^ -lmpi -L${NCDIR}/lib -I${NCDIR}/include -lnetcdff -o $@ 

pnetcdf: wam_mpi_module.o wam_file_module.o wam_timopt_module.o \
	wam_model_module.o wam_source_module.o wam_fre_dir_module.o wam_interface_module.o \
	wam_grid_module.o wam_current_module.o wam_special_module.o wam_nest_module.o \
	wam_ice_module.o wam_output_module.o wam_print_module.o wam_general_module.o \
	wam_output_set_up_module.o wam_mpi_comp_module.o wam_netcdf_metno_module.o wam_coordinate_module.o \
	read_current_input.o read_ice_input_metno.o wam_topo_module.o read_topo_input.o jafu.o \
	make_netcdf.o 
	$(FC) $(INC) $(FFLAGS) $^ -lmpi -L${NCDIR}/lib -lnetcdff -o $@ 


# dependencies:
chief.o: wam_general_module.o wam_mpi_comp_module.o
preproc.o: wam_general_module.o
wam_coordinate_module.o: wam_file_module.o
wam_general_module.o: wam_coordinate_module.o
wam_current_module.o: wam_model_module.o
wam_ice_module.o: wam_mpi_module.o
wam_output_module.o: wam_output_set_up_module.o wam_topo_module.o wam_source_module.o \
	wam_mpi_comp_module.o
wam_output_set_up_module.o: wam_special_module.o wam_timopt_module.o
wam_mpi_comp_module.o: wam_nest_module.o wam_fre_dir_module.o wam_output_set_up_module.o
wam_user_module.o: wam_assi_set_up_module.o
wam_nest_module.o: wam_grid_module.o
wam_special_module.o: wam_mpi_module.o wam_general_module.o
wamodel.o: wam_assi_module.o
wam_file_module.o: 
wam_grid_module.o:
wam_fre_dir_module.o:
wam_mpi_module.o:
wam_timopt_module.o:
wam_interface_module.o:
wam_model_module.o:
wam_topo_module.o:
wam_source_module.o:
wam_wind_module.o:
wam_boundary_module.o: wam_mpi_comp_module.o
#wam_propagation_module.o:
#preproc_module.o:
#wam_coldstart_module.o:
#wam_restart_module.o:
#wam_initial_module.o:
#wam_radiation_module.o:
#wam_assi_set_up_module.o:
read_topo_input.o: wam_general_module.o
#read_topography.o:
#read_preproc_user.o:

# rules for creating objects that need netcdf
read_wind_input_metno.o: read_wind_input_metno.f90
	$(FC) $(INC) $(FFLAGS) -I${NCDIR}/include -L${NCDIR}/lib -lnetcdff  -c -o $@ $<

read_ice_input_metno.o: read_ice_input_metno.f90
	$(FC) $(INC) $(FFLAGS) -I${NCDIR}/include -L${NCDIR}/lib -lnetcdff  -c -o $@ $<

wam_netcdf_metno_module.o: wam_netcdf_metno_module.f90
	$(FC) $(INC) $(FFLAGS) -I${NCDIR}/include -L${NCDIR}/lib -lnetcdff  -c -o $@ $<

clean:
	rm -f *.o wam *.mod pnetcdf preproc 

