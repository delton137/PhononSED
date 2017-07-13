FILES=four1 mpi dans_timer lun_management main_vars InputOutput eig_project PhononSED
OBJS=$(addsuffix .o, $(FILES))

FC= gfortran

#FFLAGS=-fpp  -O3 -C -debug -traceback
FFLAGS = -O3  -cpp

all: PhononSED.x

#compilation for debugging
debug: FFLAGS += --debug --backtrace -fbounds-check
debug: PhononSED.x

PhononSED.x : $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o  ./PhononSED.x

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o : %.F90
	$(FC) $(FFLAGS) -c $< -o $@


clean:
		rm -rf *.o
		rm -rf *.mod
