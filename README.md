# Inverse probability weighting


## Basic
Weights come from: 
+ 1. design weight (sampling)
+ 2. nonresponse weight

The weighted esitmator, given by:
$T_y $ $=\sum d_{ij} q_{ij} R_{ij} y_{ij}$

Basic parametric model (with cluster effects):
$pr(R_{ij}=1|u_i)=h(x_{ij}\beta+u_i)$
h() is the inverse logit function. 

## IPW v.s. MMI
IPW needs a model for the probability that an individual is a complete case

MI needs a model for the distribution of the missing data given the observed data

+ MI advantages:
1. Does not restrict to fully observed variables
2. Can reduce uncentrainty
3. More efficiency 

+ IPW advantages:
1. Easy to use
2. Avoid the bias from mis-specified imputation model

## A simulation 
Number of clusters in each group: 50 \\
Number of individuals in each cluster: 25\\
Repeat time: 50\\
Missing percentage: 30%
ICC: 0.06
Time Used 6.539501 hours

#### Data Generation
$Y_{ijl}=\beta_0+\beta_1 X_{ijl} +\beta_2 trt_{ijl}+ u_{ij}$
logit $R_{ijl}=\beta_0+\beta_1 X_{ijl} + u_{ij}$

i: is the indicator of treatment effect (i=1, intervention group; i=0, control group)\\
j: the jth cluster 

l: the lth individual in cluster j group i. 

#### Imputation model: 

$Y_{ijl}=\beta_0+\beta_1 X_{ijl} + u_{ij}$

#### Analysis model:
$Y_{ijl}=\beta_0+\beta_1 X_{ijl} +\beta_2 trt_{ijl}+ u_{ij}$

### Results 

Method|Package|Estimated Mean | Estiamted Std | 97% CI | Converage Rate
--- | --- | --- | --- | ---| ---
IPW|glm, lmer|4.97|0.45|(4.09,5.86)|0.96
IPW|glmer, lmer|4.99|0.45|(4.10,5.87)|0.96
Standard MI|jomo|3.38|0.57|(2.26,4.5)|0.08
MMI|jomo (cluster)|4.71|0.65|(3.44,5.99)|0.96
.|jomoImpute|4.32|0.74|(2.87,5.78)|0.94
 .|panImpute|4.34|0.74|(2.89,5.78)|0.90
 .|jomoImpute|4.99|0.62|(3.78,6.21)|0.98
 .|panImpute|5.01|0.63|(3.78,6.23)|0.98


### Discussion 
Without the consideration, the standard MI got really bad results. It underestimated the estiamtor and had a wider CI. And also under- coverage. 

Different packages for MMI has different effects and efficiency. jomo with cluster effects had the best result.  And here we can feel the advantages of IPW. 




## Two mean papers
