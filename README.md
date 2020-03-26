# Coordinated Factor Analysis (CFA)

This is an implementation of the Coordinated Factor Analysis (CFA) algorithm or also known as globally coordinated mixture of factor analysers. The main algorithm is in the function 'CFA.m'.

Three demonstrative examples are provided that can be executed via the script 'RunExamples.m'. They include a portion of the benchmark SwissRoll data set, a 2D plane of data curved in 3D with a sinusoidal signal (an S-shape), human motion capture data of a walking cycle captured with the Xsens Moven suit. An interactive figure is used to show mapping from latent variable space to the data space for each example. For the motion capture data we compare PCA and CFA reconstruction with two interactive plots. For the S-shape example there is an iterative plot to show the model as the parameters are fit.

Please refer to the following article for more details as the code is a subset of this study:

  *M. Field, D. Stirling, Z. Pan and F. Naghdy, "Learning Trajectories for Robot Programing by Demonstration Using a Coordinated Mixture of Factor Analyzers," in IEEE Transactions on Cybernetics, vol. 46, no. 3, pp. 706-717, March 2016. [doi.org/10.1109/TCYB.2015.2414277](https://doi.org/10.1109/TCYB.2015.2414277 "Learning Trajectories for Robot Programing by Demonstration Using a Coordinated Mixture of Factor Analyzers")

**Corresponding author:** Matthew Field

matthew.field@unsw.edu.au

![alt text](https://github.com/mffield/CFA/blob/master/data/CFAdemo1.PNG "Demo 1")
![alt text](https://github.com/mffield/CFA/blob/master/data/CFAdemo2.PNG "Demo 2")