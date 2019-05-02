# Engerer2-separation-model
Code for the Engerer2 diffuse fraction separation model following the newly parameterised version published in the Journal of Renewable and Sustainable Energy. The model estimates the diffuse horizontal irradiance from inputs of global horizontal irradiance, latitude, longitude and time.

The Engerer Separation Model was first described and validated in the paper: 
Engerer, N. A. 2015. Minute resolution estimates of the diffuse fraction of global irradiance for southeastern Australia. Solar Energy. 116, 215-237. 
The paper presents 3 variations of the Engerer separation model: 1, 2 & 3
The Engerer1 is fit to Non-Cloud Enhancement data, the Engerer2 to Cloud
Enhancement data and Engerer3 to clear-sky data only.

This original model was trained and tested exclusively on Australian data, and so lacked global scope, despite validating with excellent performance on many comparative studies.
This repository illustrates a new performance of the Engerer2 model as per the new paper:
Bright, Jamie M. & Engerer, Nicholas A. 2019. Engerer2: Global re-parameterisation, update and validation of an irradiance separation model at different temporal resolutions. Journal of Renewable & Sustainable Energy. 11(2), xxx.

We request that Engerer (2015), Bright & Engerer (2019), as well as this repository be cited if using this code in research. 

## PREDICTORS
The predictor variables to fit B0-5 are based on:
kt = clearness index
zen = zenith angle
ast = apparent solar time
dktc = ktc - kt (ktc is the clear-sky clearness index = ghi,cs/e_exth.
k_de = Proportion of kd attributable to cloud enhancement

## PARAMETERS
The new parameters are defined in the paper as:

 | Parameter | Predictor | Orig1-min | New1-min | 5-min | 10-min | 15-min | 30-min | 1-hr | 1-day | 
 |---|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
 | C | - | 0.042336 | 0.10562 | 0.093936 | 0.079965 | 0.065972 | 0.032675 | -0.0097539 | 0.32726 | 
 | B0 | - | -3.7912 | -4.1332 | -4.5771 | -4.8539 | -4.7211 | -4.8681 | -5.3169 | -9.4391 | 
 | B1 | K_t | 7.5479 | 8.2578 | 8.4641 | 8.4764 | 8.3294 | 8.1867 | 8.5084 | 17.113 | 
 | B2 | AST | -0.010036 | 0.010087 | 0.010012 | 0.018849 | 0.0095444 | 0.015829 | 0.013241 | 0.13752 | 
 | B3 | zen | 0.003148 | 0.00088801 | 0.003975 | 0.0051497 | 0.0053493 | 0.0059922 | 0.0074356 | -0.024099 | 
 | B4 | dKtc | -5.3146 | -4.9302 | -4.3921 | -4.1457 | -4.169 | -4.0304 | -3.0329 | 6.6257 | 
 | B5 | Kde | 1.7073 | 0.44378 | 0.39331 | 0.37466 | 0.39526 | 0.47371 | 0.56403 | 0.31419 | 

## INPUT PARAMETERS
+ ghi   = Global horizontal irradiance [W/m^2]
+ time = in the matlab datevec format [yyyy, mm, dd, HH, MM, SS] [UTC]
+ latitudes = either a latitude per measurement or one for all.
+ longitude = either a longitude per measurement or one for all.
+ averaging_period = x in minutes. Must be 1, 5, 10, 15, 30, 60 or 1440 

## OUTPUT PARAMETERS
+ Kd = diffuse fraction = Edh/ghi [fraction]
+ Edh = diffuse horizontal irradiance [W/m^2]
+ Ebn = direct/beam normal irradiance [W/m^2]
