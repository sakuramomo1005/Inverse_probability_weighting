
## Highlight


What I changed for the simulation:

* 1. Fit one simulation to compare the results from package *lme4*+*geepack* and package *CRTgeeDR* (the results are a little bit different)
* 2. For result 1, since I calculated the predicted y values by myself, I changed to use **preidct()** function this time. Use the **predict()** function to get the correct weigths from *glm* and *glmer*.
* 3. Re-run the code with *varied cluster size.* m=25 to UNI(20,30), m=50 to UNI(45,55) 
* 4. Check the boxplot of weigths


## Files

* *ipw_result2.R* contains the simulation codes and results output codes
* The *ipw_result2.pdf* is the detailed results report
* The *ipw_result2.Rmd* is the Rmarkdown file used to generate ipw_result2.pdf
* The ppt file are the meeting slides
