OBJECTS = wam_general_module.o wam_print_module.o wam_file_module.o \
wam_print_user_module.o  read_spectra_file.o \
read_spectra_user.o wam_coordinate_module.o \
wam_mpi_module.o  wam_timopt_module.o \
wam_model_module.o wam_source_module.o wam_fre_dir_module.o wam_interface_module.o \
wam_grid_module.o wam_current_module.o wam_special_module.o wam_nest_module.o \
wam_ice_module.o wam_output_module.o  \
wam_output_set_up_module.o wam_mpi_comp_module.o wam_spec_netcdf_metno_module.o  \
read_time_user.o \
read_current_input.o read_ice_input_metno.o wam_topo_module.o read_topo_input.o jafu.o \
make_spec_netcdf_metno.o

pnetcdf_spec:
        #mpif90 $(OBJECTS) -o pnetcdf_rad     \
        #-I$(NETCDFHOME)/include -L$(NETCDFHOME)/lib -lnetcdf -lnetcdff
	ifort -o pnetcdf_spec $(OBJECTS) -lnetcdff -lnetcdf -lmpi 
