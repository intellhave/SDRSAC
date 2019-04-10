#------------------------------------------------------------------------------#
#                                                                             
#                                 VARIABLES 
#                                                                             
#------------------------------------------------------------------------------#
CXXFLAGS = -O2 -Wall -fmessage-length=0
MXX = mex
MEXEXT = $(shell mexext)
ARCHIVENAME = kdtree
#NOT NECESSARY?
#MXXFLAGS = -I/Applications/matlab/extern/include

#------------------------------------------------------------------------------#
#                                                                             
#                       WHICH COMPONENTS WILL BE BUILT?
# Could also use: #TARGETS = $(wildcard cpp/*.cpp)
#------------------------------------------------------------------------------#
SOURCES = $(wildcard *.cpp)
TARGETS = $(SOURCES:%.cpp=%.$(MEXEXT))

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
all: $(TARGETS)

#--- Compiles in debug mode
debug: CXXFLAGS += -g 
debug: all

##--- Compile MEX toolbox
%.$(MEXEXT) : %.cpp $(HDRS)
	$(MXX) $(MXXFLAGS) -I/Applications/matlab.app/extern/include $<      

##--- Clean steps
clean:
	@rm -vf *.$(MEXEXT) 
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
	cp *.$(MEXEXT) $(ARCHIVENAME)
	cp *.mexmaci64 $(ARCHIVENAME)
	cp *.mexw64 $(ARCHIVENAME)
	cp *.m $(ARCHIVENAME)
	cp *.h $(ARCHIVENAME)
	cp *.cpp $(ARCHIVENAME)
	cp Makefile $(ARCHIVENAME)
	cp README $(ARCHIVENAME)
	#tar -cvf $(ARCHIVENAME).tar.gz $(ARCHIVENAME)
	zip -r -v ${ARCHIVENAME}.zip ${ARCHIVENAME}
	rm -rf $(ARCHIVENAME)
	mv ${ARCHIVENAME}.zip ../
