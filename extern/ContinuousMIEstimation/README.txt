======BASIC DESCRIPTION:=====

This software implements the modified Kraskov-Stoegbauer-Grassberger (KSG)
estimator of mutual information between two continuous variables, with
modifications introduced by Holmes and Nemenman in

1. Holmes, CM. & Nemenman, I.  (2019) Estimation of mutual information for
real-valued data with error bars and controlled bias. Phys Rev E, 100(2), 022404.

Please cite the reference [1] above when using this software. The software
is COPYRIGHTED by Holmes and Nemenman, 2019, and is distributed under
GPL 3. Please read the license information in the attached file.

The original KSG estimator was introduced in

2. Kraskov, A., Stoegbauer, H., & Grassberger, P. (2004). Estimating
mutual information. Phys Rev E, 69(6), 066138.

3. Stoegbauer, H., Kraskov, A., Astakhov, S.A., & Grassberger,
P. (2004), Least Dependent Component Analysis Based on Mutual
Information. Physical Review E, 70, 066123.

The code for the KSG estimator was copyrighted and distributed under
GPL 3 in 2009, by Astakhov, Kraskov, Stoegbauer, and Grassberger. The
original KSG code is distributed with no changes together  with this
software for ease of use and completeness. Please cite references [2]
and [3] when using this software as well.

All comments, suggestions, and proposed modifications should be
communicated to the authors:
	     Caroline M. Holmes <cholmes@princeton.edu>
	     Ilya Nemenman <ilya.nemenman@emory.edu>

The main GitHub repository for this software is:
    	     https://github.com/EmoryUniversityTheoreticalBiophysics/ContinuousEntropyEstimation
The software can be downloaded and modified there.

The software requires a C compiler with standard libraries and matlab 
compatibility, as well as a reasonably modern version of Matlab (we have 
run it using various Matlab versions R2016b-2020b, but newer and older 
versions may also work).


=====CONTENT OF THE PACKAGE=====

The contents of this folder are as follows:

1. A folder titled ‘kraskovStoegbauerGrassberger’. This contains the
files that perform the nearest-neighbor mutual information estimation
using the KSG estimator (see Refs. [2,3]). The files are distributed
with no change and copyrighted by their original authors. This code is 
NOT currently called by any of the matlab routines, and is not necessary 
in order to run the methods detailed in our 2019 paper. We distribute it 
here only for ease of adapting the code to individual's needs. Those only 
wishing to implement only our code can safely ignore or not download this 
folder.

2. In the main folder, MIxnyn_matlab.m, MIxnyn.C, and miutils.h. These are 
adapted from the code in the folder 'kraskovStoegbauerGrassberger'; the 
primary change is adding a mex function. These functions run the core nearest 
neighbor mutual information estimation.
 
In order to run these files, you must have a compatible compiler for C code 
(see https://www.mathworks.com/support/requirements/supported-compilers.html), 
and you must enter the line "mex MIxnyn.C" into matlab command line. This 
will compile all required C code.
	
3. The .m files described in the text of Ref. [1]. For examples of how
to run these files on data, see the ExampleScripts below.

	findMI_KSG_subsampling.m - This function calculates the mutual
	   information between X and Y both for the full data set and for
	   a series of nonoverlapping subsets, and outputs both that
	   information for each of the subsets. This also allows the user
	   to check for sample size dependent bias in the mutual
	  information estimate. 
	findMI_KSG_stddev.m - This function calculates the error bars
	   for the mutual information estimate, using the chi-squared
	   method described in [1], which involves extrapolating from
	   the variances at smaller N’s to the variance at the full
	   sample size. This function takes as an input the outputs of
	   findMI_KSG_subsampling.m.  
	findMI_KSG_bias_kN.m - This function calls both
	   findMI_KSG_subsampling.m and findMI_KSG_stddev.m. It
	   performs the information estimates at various values of k,
	   allowing the user to check the k-dependence of the mutual
	   information estimate, and outputs the mutual information
	   and error bars for all requested values of k. 
	reparamaterize_data.m - As discussed in Ref. [1],
	   reparamaterizing data to a gaussian can aid in mutual
	   information estimates. This function reparamaterizes the
	   input variable to a gaussian.  

3. A folder titled ‘ExampleScripts’, which contains scripts to perform
analyses similar to those in the figures of Ref. [1]. We recommend going
through these in order, as they build in complexity. The 'tutorial' branch 
contains a more open-ended informal introduction to the code, if that is 
preferred.

	a. gaussianExample.m - this function performs mutual information
	    estimation on correlated gaussian data. It generates
	    figures similar to Figs. 2, 3 in Ref. [1], which show how
	    error bars are calculated, and how mutual information
	    estimates depends on k and on N.
	
	b. logNormalExample.m - this shows how error bars are
	    calculated and how mutual information estimates depend on
	    k and N for log-normal bivariate date. It also
	    demonstrates the effects of reparameterizing the data and
	    generates the equivalent of Fig. 4 in Ref. [1].

	c. NfkappaBDataAnalysis.m -- this function performs mutual
	   information estimation on N-kappaBData.mat. It generates two
   	   plots comparing the N-dependence of mutual information
	   estimates with and without reparamaterization. The version
	   with reparamaterization is equivalent to Fig. 7 in Ref. [1].  
	
	d. NfkappaBData.mat -- this is the data file for data used in
            Fig. 7 in Ref. [1]. The data lists single cell NF-kappaB
            (P65 nuclear localization) and p-ATF-2 activation and in
            response to a doze of the TNF stimulus. See Ref. [1] and

	    	  Cheong, R., Rhee, A., Wang, C. J., Nemenman, I., &
	    	  Levchenko, A. (2011). Information transduction
	    	  capacity of noisy biochemical signaling
	    	  networks. Science 334(6054), 354–358.

	    for additional details about the data.
		  
	e.  finchDataAnalysis.m -- this function performs mutual information 
	    estimation on finchData.mat, and is equivalent to Fig. 8 in 
	    Ref. [1].

	f. finchData.mat - this file contains the data used 
	   finchDataAnalysis.m.  This data was used in Fig 8 in Ref. 
	   [1]. The data is two consecutive interspike intervals (X), 
	   and the following two consecutive interspike intervals (Y), 
	   for a motor unit in an anesthetized finch. See Ref. [1] and
            
            	  Srivastava, K., Holmes, C. M., Vellema, M., Pack, A. R., 
		  Elemans, C. P. H., Nemenman, I., Sober, S. J. (2017). 
		  Motor control by precisely timed spike patterns. Proc 
	  	  Natl Acad Sciences 114(5), 1171-1176.


=====OVERVIEW OF HOW TO USE THE PACKAGE=====


The function findMI_KSG_bias_kN.m can be used as a one-line command to 
estimate mutual information. It allows as inputs flags that determine 
whether a large number of plots are created; these plots mirror the analyses 
discussed in our 2019 paper about choosing the parameter k and checking 
for bias.

If reparameterization is required, one should reparameterize the data before 
plugging it into findMI_KSG_bias_kN.


