# Eccentric-Binary
A few python scripts for investigating the detection of eccentric binaries through gravitational waves

These scripts save images and csv files to the workspace (since Docker does not have graphics capabilities) so you may need to delete these after running the code.

chisq_final.py finds the SNR and chisq value for a selection of different chirp masses and eccentricities

newsnr.py finds the power in the n=2 mode for the different frequency bins and "fixes" the SNR by multiplying the SNR by a factor for each bin

snrtimeshift.py Plots the SNR time series for different values of eccentricity and can also investigate the addition of noise
