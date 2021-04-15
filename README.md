# multiridge
R package for multi-penalty ridge regression
NOTE: versions V1.7 and higher are upto 10 times faster than earlier versions 

library(devtools);
install_github("markvdwiel/multiridge")

You may also also install multiridge by downloading the .zip or tar.gz file, and use in R (with correct file name):
install.packages(filename, repos=NULL); install.packages(c("pROC", "risksetROC", "survival", "mgcv"))

Demo script and data available from: https://drive.google.com/open?id=1NUfeOtN8-KZ8A2HZzveG506nBwgW64e4

Manuscript: Mark A. van de Wiel, Mirrelijn van Nee, Armin Rauschenberger (2021). Fast cross-validation for high-dimensional ridge regression. J Comp Graph Stat
Use multiridge 1.3 (zip and tar.gz available from this repository) to reproduce results, as versions >=1.4 no longer depend on the package 'penalized', which was orphaned. Scripts and data used for the manuscript, plus script to check results: https://drive.google.com/drive/folders/1hwQEezOQZATZb0hixg67KXxC3Rm-2zvq?usp=sharing

