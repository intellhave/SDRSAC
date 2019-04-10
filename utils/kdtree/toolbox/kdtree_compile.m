function kdtree_compile()
% clc
localpath = fileparts(which('kdtree_compile'));
fprintf(1,'Compiling kdtree library [%s]...\n', localpath);

err = 0;
err = err | mex('-outdir',localpath, [localpath,'/kdtree_build.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/kdtree_delete.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/kdtree_k_nearest_neighbors.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/kdtree_ball_query.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/kdtree_nearest_neighbor.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/kdtree_range_query.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/kdtree_io_from_mat.cpp']);
err = err | mex('-outdir',localpath, [localpath,'/kdtree_io_to_mat.cpp']);

if err ~= 0, 
   error('compile failed!'); 
else
   fprintf(1,'\bDone!\n');
end       