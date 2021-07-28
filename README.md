# rostral-caudal-AC-ForwardModel
**Predicting Neuronal Response Properties from Hemodynamic Responses in the Auditory Cortex**

_Authors: Isma Zulfiqar*, Martin Havlicek*, Michelle Moerel, Elia Formisano_
_*equal contribution_

In the folder **Code**, you will find two sub-folders and a matlab script (.m) file.

1. **Folder Neuronal Model** contains main script _**AC_model.m**_: 

This script simulates neuronal responses using Wilson Cowan Cortical Model (WCCM, Zulfiqar et al., 2020) with varying spectro-temporal properties (indicated by tau: TIME CONSTANT and Q: QUALITY FACTOR for frequency tuning curves) in the Belt region. This is acheieved by manipulating the: 

    Time constant, tau (for temporal dynamics)
    The connectivity kernel (controlling connectivity between the simulated core and belt region), and
    Sigmas (spatial spread of activation)
    
Overall, these parameters control the dynamics of the belt region. In the script, for a single model configuration, a sample input signal (tone) is passed through the model to mimic the stimulus presented in the fMRI experiment (Santoro et al., 2017) in a fast event related design where stimulus is presented in silent interval after acquisition.

    Output of this script (for all 28 models of the Belt, configuration listed in Supplementary Table 1 of the manuscript) in response to all 12 run of randomly presented 256 sound stimuli presented in ther experiment by Santoro et al. (2017) is uploaded at '...\Hemodynamic Model\neuronal simulations\subj1'


2. **Folder Hemodynamic Model** contains 3 subfolders.

_**Masked timecourses**_ contains sample measured masked ('slow' and 'fast' belt region masks) timecourse dataset for a single participant, single hemisphere. _**stimOrder**_ contains order in which the sounds were presented to this participant, in each run (total 12 runs). _**Neuronal Simulations**_ contains output of neuronal model for all 28 models (28 tau Q pairs signifying different response properties), for all sound stimuli (12 runs, as presented by Santoro et al., 2017) simulating cortical responses for multiple belt regions (with varying tau and Q).


The main script _**Model_fitting_comparison.m**_ simulates BOLD responses from Neuronal Responses (from folder _Neuronal Simulations_) generated via WCCM with varying spectral and temporal processing in the auditory belt. The script also performs comparison with measured responses (from folder _Masked timecourses_) using VB-GLM to predict the best neuronal model for the measured response. The output is saved in a .mat file.

3.  **Output_Analysis.m** uses .mat file produced by Model_fitting_comparison.m to analyze the results and generate figure 6 (Model Predictions) and figure 9 (Distribution of voxels across simulated regions) from the main manuscript.


Contact info: isma.zulfiqar@maastrichtuniversity.nl
