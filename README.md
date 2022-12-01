## PBPK modelling of serum levels of vitamin D and 25-hydroxyvitamin D
Reference:
Huang Z & You T (2021) Personalise vitamin D3 using physiologically based pharmacokinetic modelling.
CPT: Pharmacometrics Syst Pharmacol. 10(7):723-734.
doi: 10.1002/psp4.12640.

## Significance of the work
1. The first-ever well-qualified pharmacokinetics model of vitamin D
   Data Set    Purpose     Observation    No. of Arms     Dose Range
          A   Training              VD             13     [70, 2500] μg SD; [20, 275] μg QD
          B   Training         25(OH)D             43     10 and 100 μg QD
          C       Test         25(OH)D             83     [12.5, 1250] μg QD
          D       Test         25(OH)D             16     [1250, 50000] μg SD
2. Provides accurate predictions for diverse conditions and forms a basis for personalising vitamin D dosing regimen
3. Simple structure for clear communications

## Model features
1. Accurately predicted most test data except for slight overprediction for extremely high single doses
2. Simple in structure: only considers metabolism and 
3. One main hypothesis is strongly supported: each arm’s unique endogenous vitamin D synthesis rate explains data variability
4. Reliable inferences for parameter values
5. A single set of parameter values to encompass diverse data (different countries, seasons, races, gender, age, weight, BMI etc)

## Files
1. 25(OH)Dtestset.csv  Test set: Figure S7. Visual predictive check of 25(OH)D after repeat dose.
2. 25(OH)Dpc10_100.csv Training set: Figure S6. Visual predictive check of 25(OH)D after repeat dose.
3. Final_model.R       Final model for VPC.
4. MCMC-Run009b.RDS    MCMC results needed for sampling the final model.

## To run the R script
1. Save all files into a directory on your computer
2. Set the work directory to the path where csv files are located by modifying and running line 38 of the R script
   setwd("../project/VDmodel")
3. Computation time for sampling varies by the duration you choose to simulate. Roughly speaking, 5 minutes to simulate 60 days. Have a look at Figure S6 and Figure S7 in the supplementary materials and choose which set you'd like to simulate
4. Have fun!
