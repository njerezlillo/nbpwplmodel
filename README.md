# A new regression piecewise power-law cure fraction model: diagnostic analytics and medical application

"pwpowerlaw.R": Contains all the main functions related to the piecewise power-law model.

"nbpwpl.R": Contains all the main functions related to the NBPWPL regression model.

"pwexp.R": Contains all the main functions related to the piecewise exponential model.

"nbpwexp.R": Contains all the main functions related to the NBPWEXP regression model.

"nbpwpexp.R": Contains all the main functions related to the NBPWPEXP regression model.

"nbcfwei.R": Contains all the main functions related to the cure fraction Weibull model.

## Application

"1 descriptive.R": Contains the script for all the exploratory analysis performed in the application.

"2 change_point_detection.R": Contains codes for estimating the location and optimal number of change points.

"3 nbpwpl_reg_model.R": Contains codes for estimating all parameters of the NBPWPL regression model.

"4 nbpwpl_lasso.R": Contains codes for lasso penalty and shrinkage parameter selection.

"5 nbpwpl_reg_model_final.R": Contains codes for estimating the non-penalized parameters of the NBPWPL regression model.

"6 nbpwpl_local_influence.R": Contains codes for applying local influence analysis.

"7 nbpwpl_post_deletion.R": Contains codes for the post-deletion analysis of the deleted observations.

"8 model_comparison.R": Contains codes for comparing the goodness-of-fit of all tested models.

"9 nbpwpl_mean_residual_file.R": Contains codes for computing the mean residual time function.

"kidney.csv": Corresponds to the database used.

## Simulation Study

"scenarios.R": Contains a brief exploration of the cases considered in the simulation study.

"sim_i.R": Script for generating the tables related to case i in the simulation study.

"sim_tau_i.R": Script for generating the tables related to the change points in the i-th scenario.

## Reference

Jerez-Lillo N, Tapia A, Lachos V, Ramos P. A new piecewise power-law cure fraction model: diagnostic analytics and medical application.
