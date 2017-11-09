# Inverse probability weighting
The project of applying IPW to handling missing data in CRTs

## A simulation 
Number of clusters in each group: 50 
Number of individuals in each cluster: 25
Repeat time: 50
ICC: 0.06
Time Used 6.539501 hours

#### Data Generation


#### Imputation model: 

Y_{ijl}=\beta_0+\bata_1 X_{ijl} + u_{ij}

#### Analysis model:
Y_{ijl}=\beta_0+\bata_1 X_{ijl} +\beta_2 trt_{ijl}+ u_{ij}

### Results 

Method|Package|Estimated Mean | Estiamted Std | 97% CI | Converage Rate
--- | --- | --- | --- | ---| ---
IPW|glmer, lmer|4.99|0.45|(4.10,5.87)|0.96
Standard MI|jomo|3.38|0.57|(2.26,4.5)|0.08
MMI|jomo (cluster)|4.71|0.65|(3.44,5.99)|0.96
.|jomoImpute|4.32|0.74|(2.87,5.78)|0.94
 .|panImpute|4.34|0.74|(2.89,5.78)|0.90
 .|jomoImpute|4.99|0.62|(3.78,6.21)|0.98
 .|panImpute|5.01|0.63|(3.78,6.23)|0.98


### Discussion 

Different packages for MMI has different effects and efficiency.  




## Two mean papers
