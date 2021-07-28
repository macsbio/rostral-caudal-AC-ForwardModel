# rostral-caudal-AC-ForwardModel
Predicting Neuronal Response Properties from Hemodynamic Responses in the Auditory Cortex

Authors: Isma Zulfiqar*, Martin Havlicek*, Michelle Moerel, Elia Formisano

*equal contribution

In the attached compressed zip folder, you will find two sub-folders and a script (.m) file.

1. **Folder Neuronal Model** contains main script _**AC_model.m**_: 

This script simulates neuronal responses using Wilson Cowan Cortical Model (Zulfiqar et al., 2020) with varying spectro-temporal properties (indicated by tau: TIME CONSTANT and Q: QUALITY FACTOR for frequency tuning curves) in the Belt region. This is acheieved by manipulating the: 

    Time constant, tau (for temporal dynamics)
    The connectivity kernel (controlling connectivity between the simulated core and belt region), and
    Sigmas (spatial spread of activation)
    
Overall, these parameters control the dynamics of the belt region. In the script, a sample input signal (tone) is passed throught the model to mimic the stimulus presented in the fMRI experiment (Santoro et al., 2017) in a fast event related design where stimulus is presented in silent interval after acquisition.

    Sample output of this script in response to 1 run of sound stimuli presented by Santoro et al. (2017) is located at '...\Hemodynamic Model\neuronal simulations\subj1'


2. **Folder Hemodynamic Model** contains 3 subfolders.

_**Masked timecourses**_ contains sample measured masked timecourse dataset for a single participant, single hemisphere. _**stimOrder**_ contains order in which the sounds were presented to this participant, in each run. _**Neuronal Simulations**_ contains sample output of neuronal model to 1 model (one tau Q pair), single run of sound stimuli (as presented by Santoro et al., 2017) simulating cortical responses for multiple belt regions (with varying tau and Q).


The main script _**Model_fitting_comparison.m**_ simulates BOLD responses from Neuronal Responses generated via WCCM with varying spectral and temporal processing in the auditory belt. The script also performs comparison with measured responses using VB-GLM to predict the best neuronal model for the measured response. The output is saved in a mat file.

3.  **Output_Analysis.m** uses .mat file produced by Model_fitting_comparison.m to analyze the results and generate figure 6 (Model Predictions) and figure 9 (Distribution of voxels across simulated regions) from the main manuscript.


Contact isma.zulfiqar@maastrichtuniversity.nl if you require intermediary output files.
