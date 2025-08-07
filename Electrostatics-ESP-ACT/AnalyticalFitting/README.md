AnalyticalFitting
=================
        
+ Fit_ESP-4models.py computes the electrostatic potential for each ion and determines the parameters for each charge model.

+ Electrostatic_4models.py can be used to calculate the electrostatic interactions for each model and to plot the SAPT0 electrostatics for the interactions between each ion pair.


## Limits of ESP Fitting

The use of point charges to model molecular charge distributions introduces significant errors in the electrostatic potential (ESP), as detailed in our theoretical analysis.

We show that:

- Point charges do not correctly reproduce quantum mechanical ESPs.
- ESP fitting is only valid within the fitted region and cannot generalize accurately.
- Errors in fitted charge distributions lead to significant deviations in electrostatic interaction energies.
- Even advanced models based on Slater-type charges suffer if only ESP is used for fitting.

