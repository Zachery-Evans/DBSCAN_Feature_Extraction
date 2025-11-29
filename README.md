DBScan Feature Extraction is meant to do a PCA analysis on an IR spectra formatted in Quasar. 

The program is currently set up to select 2 standard deviations from the mean of the PE peak to omit background noise, train the PCA model on that data, and then transform the total dataset (including the background noise) using the model trained
only on the selected data. 

The work still to be done includes the feature extraction from the DBScan to isolate the different latent variables the model discovers. 
