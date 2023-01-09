# glmseg
glmseg is an R segmentation package based on max-EM algorithm (an extension of the Forward Backward algorithm for constrained HMM which is estimated using EM algorithm initialized with Fused Lasso or Binary Segmentation algorithm).

In this first version, we provide the package with the homoscedastic Gaussian distribution. Other loss functions based on heteroscedastic Gaussian, Poisson, logistic regression and Weibull survival models (already functional and in the process of being documented) will follow.

To use the package, you just need to
- clone the project code in this Git repository,
- open the R project (glmseg.Rproj file) in RStudio,
- go to the "Build" menu of RStudio and click on "Install" (with the hammer symbol). 
This will load the library into RStudio, and you will be able to see more details about the input data and the output results of the glmseg_gauss_homoscd function. You will also be able to see an example of using the code with simulated data.

Please, note that the code is still under development. Do not hesitate to make constructive remarks for its improvement.
