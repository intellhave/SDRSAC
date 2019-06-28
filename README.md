# SDRSAC

MATLAB implementation of the following CVPR'19 paper:
* SDRSAC - Semidefinite-Based Randomized Approach for Robust Point Cloud Registration without Correspondences.- Huu Le, Thanh-Toan Do, Tuan Hoang, and Ngai-Man Cheung (Oral).

Paper can be accessed at: https://arxiv.org/pdf/1904.03483v1.pdf

Before running the demo, please download the latest version of SDPNAL+ (https://drive.google.com/open?id=1wE90nFu95Lq4AOazq_rTaI6D70h4NhK-) and extract the zip file  to /solvers/SDPNAL+v1.0/ (you can overwrite the existing folder)



Run ``demo_sdrsac.m `` to start the demo.

Notes: 
* Currently, SDPNAL++ (provided in the ``solver`` folder) is used as the default solver. Better solvers can be replaced to improve the performance.
* The implementation of known correspondences will be available soon.
* Tests with real-data will be added.
