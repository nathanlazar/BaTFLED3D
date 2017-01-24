BaTFLED3D
=============
Bayesian Tensor Factorization Linked to External Data in three dimensions

BaTFLED is a machine learning algorithm designed to make predictions and determine interactions in data that varies along three independent modes. For example BaTFLED was developed to predict the growth of cell lines when treated with drugs at different doses. The first mode corresponds to cell lines and incorporates predictors such as cell line genomics and growth conditions. The second mode corresponds to drugs and incorporates predictors indicating known targets and structural features. The third mode corresponds to dose and there are no dose-specific predictors (although the algorithm is capable of including predictors for the third mode if present).

BaTFLED assumes a generative model where matrices of predictors are multiplied through projection matrices to form latent representations of the input data. These latent representations vectors are analogous to principal component vectors. Matrices containing the latent vectors are combined to form the response tensor which, in this case, are measures of the growth of cell lines after treatment with different drugs at different doses.

The algorithm learns distributions on values in the projection matrices and latent matrices using a Bayesian framework. Prior distributions on the projection matrices can encourge a sparse representation of the input data. Values in the model are learned through variational approximation which maximizes a lower bound on the posterior likelihood by iteratively updating values in the model. This process is deterministic given a random initialization and update ordering.

Model diagram
-------
![](vignettes/BaTFLED_model.png)

Running BaTFLED on a simulated dataset
-------

Please see the Rmarkdown vignette 'vignettes/BaTFLED3D_vignette.html'.

Installation
-----------

```
install.packages('BaTFLED3D_0.1.5.tar.gz', repos = NULL, type="source")
library(BaTFLED3D)
```
