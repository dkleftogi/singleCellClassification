## FAQ

### Question 1: How do I resolve the dependencies ?

We provide session info file with all requested packages.

Check:

```
	mySessionConfig.sh

	and install the packages with your prefered way
```

### Question 2: How do I run the pipeline ?

We provide the following scripts:

```
implementGatingStrategy.R

hemaClassifier.R

generateFeatures.R

utilityFuncs.R
```

that implement the core parts of our workflow. 

We suggest to run the analysis in the indicated order:

```
	1. implementGatingStrategy.R : develops a gating strategy to annotate healthy cell sub populations

	2. hemaClassifier.R : use the annotated cell types to develop an LDA classifier; then use this classifier to annotate cells from leukemia patients

	3. generateFeatures.R : engineer cell type-specific features suitable for Machine Learning modelling using DREMI scores

	4. utilityFuncs.R : different utility functions that must be loaded using source(utilityFuncs.R) command
	
```

### Question 3: What the pipeline is doing and what are the outputs ?

If you execute the previous codes correctly, our pipeline generates several plots for outputs. 

The plots summarise most of the main and supplementary figures presented in our publication.

They are also useful for the user to gain full control over the data. 

The generateFeatures.R code produces as output a list of cell-type specific feature vectors. It can be used in the future for Machine Learning modelling and feature selection.

Please note that the projec tis under development and more interactive vignette will be released soon. 

For now we provide a minial example of the deployed workflow.

### Question 4: Is it possible to run specific parts of the analysis and not the pipeline from scratch ?

The short answer to this question is NO. Different functions rely on input from the previous ones.

### Question 5: I have problems to run the analysis, what should I do ?

Please read carefully the error messages and try to resolve any problem with dependencies. There might be also some problem with the R studio versions, or operating system specific problems. Note that the analysis have not been tested under Windows or Linux operating systems.  

You can always contact Dimitrios and ask for help, more info and datasets are available uppon reasonable request =)


