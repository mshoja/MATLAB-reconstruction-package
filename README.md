This repository contains a MATLAB package called the MATLAB Reconstruction Package. This package implements a maximum-likelihood estimation procedure to fit the following stochastic Langevin systems to univariate time series data with different resolutions (this process is known as 'system reconstruction'):

                                    dx = mu(x) dt + sigma(x) dW,
                                                                                                                                 
where mu(x) is the drift function, sigma(x) is the diffusion function, and 𝑑𝑊 represents Brownian noise, which is uncorrelated and Gaussian. The package can reconstruct both high-resolution and low-resolution data. For high-resolution data, the Euler method (termed 'Euler reconstruction') is implemented,
while for low-resolution data, a more accurate methodology pioneered by Aït-Sahalia is implemented (termed 'Hermite reconstruction'). For details of both methods, please refer to the publication: Aït‐Sahalia, Y. (2002). Maximum likelihood estimation of discretely sampled diffusions: a closed‐form approximation approach. Econometrica, 70(1), 223-262.

This package supports two different reconstruction procedures: parametric reconstruction and spline reconstruction (for both low and high-resolution data).

Parametric Reconstruction: The user specifies a parametric model for the drift and diffusion functions (e.g., mu(x) = a*x^3-b*x^2, sigma(x) = c) and the package estimates the model parameters.
Spline Reconstruction: The user does not need to specify any model. Instead, a coarse mesh across the state space, called a 'knots sequence', should be specified. Spline modeling is more accurate, especially for Hermite reconstruction which
is used for low-resolution data, and is useful when choosing a model is challenging. 

This package includes:

A paper titled 'Reconstructing Langevin Systems from High and Low-Resolution Time Series using Euler and Hermite Reconstructions'
An appendix with solid mathematical descriptions
A user-friendly tutorial

Recommended Learning Strategy
Depending on your needs, follow these steps to learn how to use the package:

Quick Analysis:
Briefly review the paper, which contains five case studies (three simulated data sets and two real data sets from ecology and climate science). Execute the multi-section code 'AllFigures.m' to reproduce the examples in the paper. Then, use the codes for your own data.

Intermediate Learning:
Run the codes in the folder called 'Examples', which contains 10 case studies for various datasets.

In-Depth Learning:
Read the tutorial, which is explained in a user-friendly manner and includes 16 case studies for a range of simulated and real datasets.

Advanced Learning:
Study the appendix carefully for a more in-depth understanding of the mathematics, and refer to the references provided.

Usage and Redistribution
Use and redistribution of the codes in this repository are allowed, provided that you cite the following paper:

Babak M. S. Arani, Stephen R. Carpenter, Egbert H. van Nes, Reconstructing Langevin Systems from High and Low-Resolution Time Series using Euler and Hermite Reconstructions.
