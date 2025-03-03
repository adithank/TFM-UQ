# TFM-UQ
MATLAB code for 'Uncertainty-Aware Traction Force Microscopy'.  https://www.biorxiv.org/content/10.1101/2024.07.05.602172v2.abstract

## 1) PIV-UQ demo : 'demo_PIVUQ.m'

Please setup the PIV parameters in 'PIVConfig.json' with the following fields :

```json
{
  "Performance": {
    "MaxRAM": 16, % Max RAM to be used for array creation
    "UseGPU": true, % Use GPU if available
    "MaxGPUMem": 6, % Max GPU memory 
    "Precision": "single" % single or double precision variables
  },
  "Deformation": [
    {
      "wdw_size": 64, % PIV interrogation window size
      "wdw_spacing": 32, % PIV interrogation window spacing 
      "Nreps": 25, % PIV-UQ bootstrap repetitions
      "SubPixelInterpolationMode": "gaussian" % Gaussian or polynomial sub-pixel interpolation scheme
    }
    
  ]
}

```

PIV.output is a struct variable that contains the following fields : 
```matlab

xvec,yvec % X and Y coordinates of sizes (nY,nX) respectively
UPost, VPost % Post-processed PIV vector (Use this)
U, V % Unprocessed PIV vector 
Usamp, Vsamp % Bootstrap samples with a size (nY,nX,nB). nB is the number of bootstrap reps.
UStd, VStd % Standard deviation calculated from bootstrap reps.
delInds % "Bad windows" as a mask

```

## 2) Run 'MCMCGibbs_TFMUQ.m' for TFM-UQ analysis using the PIV-UQ data

