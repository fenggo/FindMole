# Makefile for GNU Linux / Debian stubs

SHELL = /bin/sh

# System-specific settings

FC=gfortran
#FC=ifort
CFLAG= -c

# Link target
FOR= modules.f90  \
     setpara.f90 \
     initial_uspl.f90 \
     info.f90 neighbor.f90 \
     findmole.f90  \
     string.f90 read.f90 write.f90 \
     pkenergy.f90 \
     memory.f90 \
     main.f90 \
     find_reac.f90 \
     component.f90 \
     product.f90 \
     distribution.f90 \
     msd.f90 lattice.f90
    

OBJ= $(FOR:.f90=.o)

# Link rule


FindMole: $(OBJ)
	$(FC) $(OBJ) -static -o findmole


# Compilation rules
%.o: %.f90
	$(FC) $(CFLAG) $<


.PHONY : clean
clean:
	rm -f *.o *.mod

    master 

