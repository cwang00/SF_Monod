objects = \
cbin_write.o  mat_dif_init.o  rand_mod.o  slim_fast.o \
slim_richards.o  min_react_lu.o  v_calc.o vtk_write.o \
NTransport.o Particles.o svode.o gnuplot_write.o

sources = \
src/cbin_write.f90  src/mat_dif_init.f90  src/rand_mod.f90 \
src/NTransport.f90 src/Particles.f90 \
src/slim_fast.f90 \
src/slim_richards.f90  src/min_react_lu.f90  src/v_calc.f90 src/vtk_write.f90 \
src/svode.f src/gnuplot_write.f90



.SUFFIXES:
.SUFFIXES: .F90 .o .f90 .f

F90 = pgf90
#F90 = gfortran

SLIM_r2.0.exe : $(source)
	$(F90) -g -o bin/SLIM.exe $(sources)
#	$(F90) -O3 -o bin/SLIM.exe $(sources)

.o : 
	$(F90) -c src/$(objects)
#.f90.o :
#	gfortran -c src/$(objects)

.PHONY : clean
clean : 
	rm -f bin/SLIM.exe $(objects)

debug :
	$(F90) -gDD -o bin/SLIM.exe $(objects)

