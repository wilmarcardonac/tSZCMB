FC     	= gfortran
LC     	= $(FC)
EXE    	= computebMz
DIR	= /home/cardonac/projects/project-tSZ-CMB-lensing-bispectrum/Code
F_FL   	= -O3 -I$(DIR) -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -g

#####################
OBJ = $(DIR)/fiducial.o $(DIR)/arrays.o $(DIR)/functions.o computebMz.o 

def:	$(OBJ)
	$(LC) $(F_FL) $(OBJ) -o $(EXE)

%.o:	%.f90
	$(FC) $(F_FL) -c $<

%.o:	%.F90
	$(FC) $(F_FL) -c $<

%.o:	%.f
	$(FC) $(F_FL) -c $<

clean :
	rm -f *.o *.mod *.ini *~  fort.* *.out *.il $(EXE)

### put dependencies here ###

computebMz.o :	$(DIR)/fiducial.o $(DIR)/arrays.o $(DIR)/functions.o 
