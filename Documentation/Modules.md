k# Modules

NIT is made of a set of different modules:
 - **II**: module for calculation of Intersection Information and Partial Information Decomposition (based on: *A. Makkeh, D.O. Theis, R. Vicente, BROJA-2PID: A cone programming based Partial Information Decomposition estimator, Entropy 2018, 20(4), 271-291; [doi:10.3390/e20040271](http://dx.doi.org/10.3390/e20040271).*  https://github.com/Abzinger/BROJA_2PID)
 - **MI**: for mutual information calculation
    - **BinnedMethods**: module for calculation of MI (and related quantities) using binned methods (based on: *Magri C, Whittingstall K, Singh V, Logothetis NK, Panzeri S: A toolbox for the fast information analysis of multiple-site LFP, EEG and spike train recordings. BMC Neuroscience 2009 10(1):81;*)
    - **ContinuousMethods**: modules for calculation of MI using copulas:
       - **MVC**: Mixed Vine Copulas for single and multi-variate MI calculation (*A. Onken and S. Panzeri (2016).  Mixed vine copulas as joint models of spike counts and local field potentials.  In D. D. Lee, M. Sugiyama, U. V. Luxburg, I. Guyon and R. Garnett, editors, Advances in Neural Information Processing Systems 29 (NIPS 2016), pages 1325-1333.* https://github.com/asnelt/MixedVineToolbox)
       - **NPC**: Non Parametric Copulas for single-variage MI calculation (*Safaai, H., Onken, A., Harvey, C. D., & Panzeri, S. (2018). Information estimation using nonparametric copulas. Physical Review E, 98(5), 053302.* https://github.com/houman1359/NPC_Info)
 - **SML**: contains wrappers for external supervised machine learning tools:
    - **SVM**: Support Vector Machine MATLAB interface for *libsvm* (*C.-C. Chang and C.-J. Lin. LIBSVM : a library for support vector machines. ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011.* https://www.csie.ntu.edu.tw/~cjlin/libsvm/)
    - **GLM**: interface for Generalized Linear Model library  *glmnet* (*Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, R. and Simon, N. http://www.stanford.edu/~hastie/glmnet_matlab/*)
 - **UML**: contains unsupervised machine learning tools:
 	- **PCA**: tools for Principal Component Analysis
 	- **kMeans**: tools for k-Means clustering
 - **Tools**: contains accessory functions
