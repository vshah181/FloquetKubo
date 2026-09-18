include make.inc
PROGRAM := floquet_kubo

#-----------------------------------MODULES------------------------------------
MOD_BZUT := bz_utilities.o
MOD_ARRY := array_utilities.o
MOD_MPIU := mpi_utilities.o
MOD_PARM := parameters.o
MOD_READ := read_files.o
MOD_WRTE := write_files.o
MOD_CNST := constants.o
MOD_KIND := kinds.o
MOD_HMLT := hamiltonian.o
MOD_TRAN := transport.o
MOD_LINT := lapack_interfaces.o
MOD_FLOQ := floquet.o
MOD_STAT := statistical_distributions.o

#-----------------THE FOLLOWING MODULES HAVE NO DEPENDENCIES-------------------
SOL_MODS := $(MOD_KIND) $(MOD_LINT)
#-------------------THE ABOVE MODULES HAVE NO DEPENDENCIES---------------------
ALL_MODS := $(SOL_MODS) $(MOD_BZUT) $(MOD_ARRY)\
	    $(MOD_PARM) $(MOD_READ) $(MOD_WRTE)\
	    $(MOD_CNST) $(MOD_HMLT) $(MOD_TRAN)\
		$(MOD_FLOQ) $(MOD_MPIU) $(MOD_STAT)


#-------------------------------------MAIN-------------------------------------
OBJ_MAIN := main.o

#----------------------------COMPILE THE PROGRAMME-----------------------------
all: $(PROGRAM)

$(PROGRAM): $(ALL_MODS) $(OBJ_MAIN) 
	$(FF) $(FFLAGS) $^ $(FLIB) -o $@

#---------------------------PREPARE THE MAIN OBJECT----------------------------
$(OBJ_MAIN): %.o: %.f90 $(ALL_MODS)
	$(FF) -c $(FFLAGS) $(FLIB) $<

#------------------------------PREPARE THE MODULES-----------------------------
$(SOL_MODS):  %.o: %.f90
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_MPIU): %.o: %.f90 $(MOD_KIND)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_BZUT): %.o: %.f90 $(MOD_KIND) $(MOD_PARM) $(MOD_ARRY)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_ARRY): %.o: %.f90 $(MOD_CNST) $(MOD_KIND)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_PARM):  %.o: %.f90 $(MOD_READ) $(MOD_ARRY) $(MOD_KIND)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_READ):  %.o: %.f90 $(MOD_CNST) $(MOD_KIND)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_WRTE):  %.o: %.f90 $(MOD_PARM) $(MOD_KIND)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_CNST):  %.o: %.f90 $(MOD_KIND)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_HMLT):  %.o: %.f90 $(MOD_PARM) $(MOD_KIND) $(MOD_CNST)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_TRAN):  %.o: %.f90 $(MOD_PARM) $(MOD_KIND) $(MOD_HMLT) $(MOD_ARRY) $(MOD_CNST) $(MOD_LINT) $(MOD_STAT)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_STAT):  %.o: %.f90 $(MOD_KIND) $(MOD_CNST)
	$(FF) -c $(FFLAGS) $(FLIB) $<

$(MOD_FLOQ):  %.o: %.f90 $(MOD_PARM) $(MOD_KIND) $(MOD_CNST)
	$(FF) -c $(FFLAGS) $(FLIB) $<

#------------------------------------------------------------------------------
clean:
	rm *.mod *.o $(PROGRAM)
