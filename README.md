# Applying MCT to disease ecology.

Instructions for code review:

The code for our decompositions are broken into four primary R scripts: "decompositions.R", "between_host_model.R", "within_host_model.R", and "model_functions.R". The script "rmpi.R" is the code to run the decompositions on Teton, whereas "parameters.R" and "analysis.R" hold the different parameter specfications depending on the scenario of interest (e.g. sensitivity analysis, neutral, and comp/col) and the code to produce the figures, respectively. You can ignore "rmpi.R" and "analysis.R" for the time being (the latter has not been cleaned yet).

"decomposition.R" can be thought of as the main script that calls the other model components (between and within host dynamics) and model functions. I would suggest working your you way through "decomposition.R" and then reviewing the other scripts as you encounter new functions (although "within_host_model.R" is straightforward enough to where you could review that before anything else). The functions in "model_functions.R", which includes the functions to calculate GRWR, should be listed in the order you encounter them. If you have the time and patience, I think it would be worth it for you to run the code and examine the R objects that are called in the different functions to ensure that nothing unexpected is occuring. I have adjusted the model parameters such that the simulations shouldn't take forever to complete.

Please let me know if you have questions!





