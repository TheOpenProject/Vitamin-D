## Vitamin-D
PBPK Modelling to predict serum levels of vitamin D and its active metabolite 25-hydroxyvitamin D.

Reference:
Huang Z & You T (2021) Personalise vitamin D3 using physiologically based pharmacokinetic modelling.
CPT: Pharmacometrics Syst Pharmacol. 10(7):723-734.
doi: 10.1002/psp4.12640.

# Files:
1. 25(OH)Dtestset.csv  Test set: Figure S7. Visual predictive check of 25(OH)D after repeat dose.
2. 25(OH)Dpc10_100.csv Training set: Figure S6. Visual predictive check of 25(OH)D after repeat dose.
3. Final_model.R       Final model for VPC.
4. MCMC-Run009b.RDS    MCMC results needed for sampling the final model.

To run the R script, you need to 

1) Save all files into a directory on your computer
2) Set the work directory to the path where csv files are located by modifying and running line 38 of the R script
   setwd("../project/VDmodel")
3) Computation time for sampling varies by the duration you choose to simulate. Roughly speaking, 5 minutes to simulate 60 days. Have a look at Figure S6 and Figure S7 in the supplementary materials and choose which set you'd like to simulate
4) Have fun!
