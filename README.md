# MR Thermometry algorithm: Near Harmonic 2D Reconstruction

This code repository implements an algorithm for MR thermometry in aqueous tissues.
The algorithm provided in this repository is based on the near-harmonic 2D reconstruction (proposed by Salomir et al. [1]) that determines the drift field in fatty regions and subsequently grows that initial estimation to the entire field of view using near-harmonic 2D reconstruction.

***
The following files are relevant for utilization of the algorithm:
* __RunSalomir.m__: Example script to run algorithm on simulated data
* __Salomir_fit.m__: Main function script
***
[1] Salomir, R., Viallon, M., Kickhefel, A., Roland, J., Morel, D. R., Petrusca, L., ... & Gross, P. (2011). Reference-free PRFS MR-thermometry using near-harmonic 2-D reconstruction of the background phase. IEEE transactions on medical imaging, 31(2), 287-301.
