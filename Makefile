# fortran compiler = intel fortran
fcomp = ifort
#switch = -module obj -w -f77rtl -fpp -D_BATCH_ -D_NOLIC_ -D_LINUX_ -D__INTEL_ -D__SHMDEBUG_ -O3 -qopenmp -parallel -intconstant -zero -assume:byterecl -check:pointers -check:bounds -check:uninit

#switch = -module obj -w -f77rtl -fpp -D_BATCH_ -D_NOLIC_ -D_LINUX_ -D__INTEL_ -D__SHMDEBUG_ -O3 -qopenmp -parallel -intconstant -zero -assume:byterecl -check:pointers -check:bounds -check:uninit


#FOR DEBUGGING
switch = -module  obj -w -f77rtl -standard-semantics -fpp -D_BATCH_ -D_NOLIC_ -D_LINUX_ -D__INTEL_ -D__SHMDEBUG_ -O0  -init=zero -init=arrays -intconstant  -assume:byterecl -heap-arrays -g -check all -check noarg_temp_created -fpe0 -traceback -debug extended 

baseList = main.for
baseObj = $(baseList:.for=.o)

modList = module.for helpers.for
modList2 = Point_Module.f90  Vector_Module.f90 Geom_Module.f90 \
SortSearch_module.f90 Qtree_module.f90  Qtree.f90
modObj = $(modList:.for=.o)
modObj2 = $(modList2:.f90=.o)

OBJPATH= ./obj

all: dir modules modules2 base Quadtree
	
dir:
	@mkdir -p $(OBJPATH)

base: $(baseObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@) -mkl

ele: $(eleObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@) -mkl

mat: $(matObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@) -mkl

modules: $(modObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@) -mkl

modules2: $(modObj2)
%.o: %.f90
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@) -mkl


rwth: $(rwthObj)
%.o: %.for
	$(fcomp) $(switch) -c $< -o obj/$(notdir $@) -mkl

Quadtree: modules modules2 base
	$(fcomp) -mkl $(switch) obj/*.o -o Quadtree
clean:
	@echo "target clean: delete obj/*"
	@rm -rf $(OBJPATH)
	
