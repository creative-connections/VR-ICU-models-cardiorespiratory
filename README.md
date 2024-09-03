# VR-ICU-models-cardiorespiratory
Cardiorespiratory models in Modelica to support simulation in VR ICU and other.

## models of respiration, mechanical ventilation and blood gas exchange
Choose normal model files to load dependent libraries first and all 3 models. Choose total model files if you want to open just one model with all dependencies included.


## normal model files

Open the dependent packages first
```
Chemical.mo
Physiolibrary.mo
```
And then open the model package
```
modelECMORespiratoryVR.mo
```

## total model files

In order to open just one model and all it's dependencies, open the selected total model file, which contains all necessary dependencies

```
total/BloodyMaryPPG2Total.mo
total/HemodynamicsRegulatedHRTotal.mo
total/LungVentilatorSCMV2Total.mo
```
