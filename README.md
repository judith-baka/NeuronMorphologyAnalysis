# NeuronMorphologyAnalysis
Script package to analyze and plot neuronal reconstructions coming from the Janelia Mouselight pipeline. The script heavily lean on the SNT FIJI (ImageJ) toolbox.
## Install
```
conda create -n nma python=3.8 pyimagej openjdk=8 -c conda-forge
conda activate nma
conda install numpy matplotlib scipy pynrrd h5py google-api-python-client gspread oauth2client umap-learn -c conda-forge
pip install rastermap miniball
pip install -e .
```
