# 20161214_hydropatterning
R scripts for estimation of growth-sustained tissue water potentials, and data sets used to train a predictive model relating potentials to lateral root patterning

Note that an interactive version of the model can be accessed as a Shiny app here:
https://nrobbins.shinyapps.io/20161105_hydropatterning_app/


File Descriptions


/hypocotyl_water_potential_calculator and /root_water_potential_calculator contain files related to estimation of tissue water potentials for growing hypocotyls and roots.

-coefficients_solutions_and_matrices.xlsx: Summary of matrices used to solve the systems of equations for calculating water flow rates and water potentials.

-potentials_output.csv: Example output from water potential calculator script.

-REGR.csv: Example relative elemental growth curve used as an input for the script.
-water_potential_calculator.R: Script for calculating tissue water potentials.

/regression_model_data_sets contains files used to build the predictive model relating tissue water potentials to lateral root patterning.

-all_growth_data.csv: Full table of relative elemental growth curves for every seedling examined as part of the manuscript. Last 4 columns are the REGR values for each of the 4 biological replicates for each experimental condition.

-regression_model_Oh7B_data.csv: Table of water potentials and lateral root patterning data from comparisons of the two maize inbreds B73 and Oh7B.

-regression_model_training_data.csv: Table of water potentials and lateral root patterning data used to train the regression model. These are from seedlings treated with varying concentrations of mannitol.

-regression_model_validation_data.csv: Table of water potentials and lateral root patterning data used to validate the regression model. These are from seedlings treated with different growth inhibitors (sodium orthovanadate, DES, citric acid, and low temperature).
