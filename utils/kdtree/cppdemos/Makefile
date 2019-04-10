#------------------------------------------------------------------------------#
#                                                                             
#                                 VARIABLES 
#                                                                             
#------------------------------------------------------------------------------#
CXXFLAGS = -O2 -Wall -fmessage-length=0
MXX = mex
MEXEXT = $(shell mexext)
MXXFLAGS = -I/Applications/matlab/extern/include
ARCHIVENAME = kdtree
EXEEXT = exe

#------------------------------------------------------------------------------#
#                                                                             
#                       WHICH COMPONENTS WILL BE BUILT?
# Could also use: #TARGETS = $(wildcard cpp/*.cpp)
#------------------------------------------------------------------------------#
TARGET += kdtree_build
TARGET += kdtree_delete

TARGET += kdtree_nearest_neighbor
TARGET += kdtree_k_nearest_neighbors
TARGET += kdtree_ball_query
TARGET += kdtree_range_query

TARGET += kdtree_print
TARGET += kdtree_load
TARGET += kdtree_save

#EXETARGET = $(TARGET:%=%.$(EXEEXT))
MEXTARGET = $(TARGET:%=%.$(MEXEXT))

#------------------------------------------------------------------------------#
#                                                                             
#                              DEPENDENCY RULES                                
#                                                                             
# Makefiles rules 
# anything (%) that terminates in .bin for which an explicit rule is not        
# available is made dependable on the file which has same name but .cpp       
# extension. The compiler (CXX) with options (CXXFLAGS) is called on each of  
# the elements that trigger the rule ($@, which is left side of ":") and      
# produces an output with filename expressed by the "first" of elements from  
# which it depends ($< or right side of ":")                                  
#------------------------------------------------------------------------------#
HDRS = KDTree.h MyHeaps.h

#--- Default rule (called when you just "make")
mex: $(MEXTARGET)
all: $(EXETARGET) $(MEXTARGET)

#--- Compiles in debug mode
debug: CXXFLAGS += -g -D DEBUG 
debug: all

#--- Compile CPP examples
%.$(EXEEXT) : cpp/%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -I./ -o ../build/$@ $<

##--- Compile MEX toolbox
%.$(MEXEXT) : %.cpp $(HDRS)
	$(MXX) $(MXXFLAGS) -I./ -o $@ $<      

##--- Clean steps
clean:
	@rm -vf ../toolbox/*.$(MEXEXT) 
	@rm -vf ../build/*.$(EXEEXT)       
	@rm -vrf *.dSYM
	@echo "---------- CLEAN COMPLETED ---------"
	
#------------------------------------------------------------------------------#
#                                                                             
#                   CREATE A DISTRIBUTION IN A ZIP FILE                       
#                                                                             
# move resources to a kdtree folder, tar them, then remove the temp directory  
# and its content completely                                                  
#------------------------------------------------------------------------------#
dist:   
	rm -rf $(ARCHIVENAME)
	mkdir $(ARCHIVENAME)
	cp *.mexmaci $(ARCHIVENAME)
	cp *.m $(ARCHIVENAME)
	cp -L *.h $(ARCHIVENAME)
	cp *.cpp $(ARCHIVENAME)
	cp CHANGES $(ARCHIVENAME)
	cp TODO $(ARCHIVENAME)
	cp README $(ARCHIVENAME)
	cp Makefile $(ARCHIVENAME)
	#tar -cvf $(ARCHIVENAME).tar.gz $(ARCHIVENAME)
	zip -r -v ${ARCHIVENAME}.zip ${ARCHIVENAME}
	rm -rf $(ARCHIVENAME)
	echo $(VAR)
