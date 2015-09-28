FC     	= gfortran
LC     	= $(FC)
EXE    	= tSZ
FITSDIR = #/home/wilmar/usr/local/lib/lib
LIBFITS = #cfitsio
INCDIR	= ./Libraries/fgsl-1.0.0
IDIR	= ./Libraries/fgsl-1.0.0
LIBDIR	= ./Libraries/fgsl-1.0.0
LDIR	= #/home/wilmar/additional-software/Healpix_3.00/lib
F_FL   	= -O3 -I$(INCDIR) -I$(IDIR) -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -g
LIB_FL 	= -L$(LIBDIR) -lfgsl #-llapack -lblas #-L$(LDIR) -L$(FITSDIR) -lhealpix -lchealpix -lhpxgif -l$(LIBFITS) -Wl,-R$(FITSDIR)
#####################
OBJ   = arrays.o fiducial.o functions.o tSZ.o 

def:	$(OBJ) $(OBJNR) $(OBJODE)
	$(LC) $(F_FL) $(OBJ) $(OBJNR) $(OBJODE) -o $(EXE)  $(LIB_FL)

%.o:	%.f90
	$(FC) $(F_FL) -c $<

%.o:	%.F90
	$(FC) $(F_FL) -c $<

%.o:	%.f
	$(FC) $(F_FL) -c $<

clean :
	rm -f *.o *.mod *.ini *~  fort.* *.out *.il $(EXE)

### put dependencies here ###

tSZ.o :	arrays.o fiducial.o functions.o 
