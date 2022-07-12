Code to accompany Stochastic Average Models paper, currently available at (arXiv preprint URL here)

### DISCLAIMER: This is research code, and comes with absolutely no guarantees. Please email me directly at mmenickelly@anl.gov for help/questions. 

### IMPORTANT CREDITS: The top level directory of this code contains minq5. This code belongs to Arnold Neumaier and can be downloaded here, but be aware version 5 is no longer maintained:

https://arnold-neumaier.at/software/minq/ 

We only include it here because we have made some edits for our own purposes, and this is the easiest way to distribute those edits. 

The top level directory also contains various subroutines of POUNDERS, which will be publicly available soon. However, formquad_indep.m is non-standard in POUNDERS, and is sufficiently different, meaning it cannot just be an external dependcy. POUNDERS belongs to Stefan M. Wild, a coauthor on this paper. 

To simply see a sample run of the code and get a sense of how to call SAM_POUNDERS on your own problems, use the just_one_run.m function with the following snippet:test_function = 'rosenbrock'; % Uses the Rosenbrock test function
data_type = 'imbalanced'; % Uses imbalanced mode of data generation
macro_seed = 1; micro_seed = 1; % Sets a random seed for problem generation
b = 1; % Specifies that the resource size is 1
num_epochs = 100; % allow for 100 effective passes through the component functions
m = 16; % For this problem, specifies that problem dimension and number of component functions is 16

Instructions to run all of the experiments illustrated in the paper:

First, for the SAG experiments, you obvious need to install SAG (as suggested, we used a mex of the C code).
This is not our method, and belongs to Mark Schmidt. It can be downloaded here:
https://www.cs.ubc.ca/~schmidtm/Software/SAG.html

run_experimentX.m 

functions, where the X is integers 1:6.
For the paper, we ran all integer macro_seeds 1:30. 

This will generate a lot of mat files and put them into a results/ directory within the tests/ directory. 

You should then run

paper_figures(1:13)

to generate all the figures in the paper. 



 
