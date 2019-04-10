# kdtree
(Automatically exported from code.google.com/p/kdtree-matlab)

This library provides a minimalist implementation of a kd-tree data structure. The implementation can be used either inside MATLAB by means of MEX calls, or as a standalone tool, directly using its C/C++ interface. The image on the side has been created within matlab from "fulltest.m". 

<div align="center">
    <img src="/media/teaser.png" width="50%"/>
</div>

## Basic API

  -   kdtree_build: k-d tree construction O( n log^2(n) )
  -   kdtree_delete: frees memory allocated by kdtree 
  -   kdtree_nearest_neighbor: nearest neighbor query (for one or more points)
  -   kdtree_k_nearest_neighbors: kNN for a single query point
  -   kdtree_range_query: rectangular range query
  -   kdtree_ball_query: queries samples within distance delta from a point

## File Structure
Everyone of the scripts/functions is complete of the following: 

 -  *.cpp: the mex implementation of the sources 
 -  *.mexmaci: the compiled version of the mex (intel mac) 
 -  *.m: the comments that you can browse with the "help" command 
 -  *_demo.m: demo file to illustrate the behavior 

## Matlab Class Wrapper
I am in the process of building a MATLAB class wrapper for the library. These are the advantages:

 - the memory will be freed automatically
 -   ability to load/save kdtree structures to/from a .mat file
 -   better organization
 
# How to compile/install in Matlab
Compiling can be done in two ways. The first is directly inside MATLAB. You can compile manually each of the files by calling the command mex within the kdtree folder from the MATLAB command line. For example: 

    >> mex kdtree_build.cpp 

Alternatively, if you are in a unix environment, you might also be able to use the provided makefile. In order to do this you need to change some of the environment variables in order to make them point to your local MATLAB installation. 

NOTE: Correctly setup I assumed you have a correctly configured MEX environment.

# Hacking and development

The *.cpp files contain a MEX interface for MATLAB. At the same time, a rich set of examples which run as standalone C++ applications is provided.  In order to compile them independetly from the MEX environment, a preprocessor condition -D CPPONLY is required. The makefile uses this flags and compiles sources in both environments: C++ and MEX. 

