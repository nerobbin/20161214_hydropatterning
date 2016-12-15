# 20161214_hydropatterning
R scripts for estimation of growth-sustained tissue water potentials, and data sets used to train a predictive model relating potentials to lateral root patterning

Note that an interactive version of the model can be accessed as a Shiny app here:
https://nrobbins.shinyapps.io/20161105_hydropatterning_app/


File Descriptions


/hypocotyl_water_potential_calculator and /root_water_potential_calculator

These folders contain files related to estimation of tissue water potentials for growing hypocotyls and roots.

-coefficients_solutions_and_matrices.xlsx: Summary of matrices used to solve the systems of equations for calculating water flow rates and water potentials.

-potentials_output.csv: Example output from water potential calculator script.

-REGR.csv: Example relative elemental growth curve used as an input for the script.
-water_potential_calculator.R: Script for calculating tissue water potentials.

/regression_model_data_sets

This folder contains files used to build the predictive model relating tissue water potentials to lateral root patterning.

-all_potentials_and_growth_curves.csv: Full table of relative elemental growth curves and estimated tissue water potentials for every seedling in the training and validation data sets.

-regression_model_training_data.csv: Table of water potentials and lateral root patterning data used to train the regression model. These are from seedlings treated with varying concentrations of mannitol.

-regression_model_validation_data.csv: Table of water potentials and lateral root patterning data used to validate the regression model. These are from seedlings treated with different growth inhibitors (sodium orthovanadate, DES, citric acid, and low temperature).
