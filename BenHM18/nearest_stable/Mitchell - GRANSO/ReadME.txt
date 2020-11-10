GRANSO is a recent solver for non-smooth optimization; see below for more information. It is more effective than HANSO that we used in our paper for comparison. 

You can run "run_test" to test GRANSO on a small nearest stable matrix problem. (Note that for large input matrices, our fast gradient algorithm works better.)  

We thank Tim Mitchell for providing us with this code. 

GRANSO can be downloaded from 
http://www.timmitchell.com/software/GRANSO/





GRANSO stands for GRadient-based Algorithm for Non-Smooth Optimization.  If you want to learn more about the underlying method, see this paper:

@article{doi:10.1080/10556788.2016.1208749,
author = {Frank E. Curtis and Tim Mitchell and Michael L. Overton},
title = {A BFGS-SQP method for nonsmooth, nonconvex, constrained optimization and its evaluation using relative minimization profiles},
journal = {Optimization Methods and Software},
volume = {32},
number = {1},
pages = {148-181},
year = {2017},
doi = {10.1080/10556788.2016.1208749},
URL = {
http://dx.doi.org/10.1080/10556788.2016.1208749
},
eprint = {
http://dx.doi.org/10.1080/10556788.2016.1208749
}
}

