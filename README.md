# SDRSAC

MATLAB implementation of the following CVPR'19 paper:
* SDRSAC - Semidefinite-Based Randomized Approach for Robust Point Cloud Registration without Correspondences.- Huu Le, Thanh-Toan Do, Tuan Hoang, and Ngai-Man Cheung (Oral).

Paper can be accessed at: https://arxiv.org/pdf/1904.03483v1.pdf



Run ``demo_sdrsac.m `` to start the demo.

Notes: 
* Currently, SDPNAL++ (provided in the ``solver`` folder) is used as the default solver. Better solvers can be replaced to improve the performance.
* The implementation of known correspondences will be available soon.
* Tests with real-data will be added.
