This code implements the algorithms to find the nearest stable matrix to a given unstable one that are described in 

On computing the distance to stability for matrices using linear dissipative Hamiltonian systems, N. Gillis and P. Sharma. 

You can run the file RunMe.m to have a comparison of these algorithms on a small 10-10 Grcar matrix. 

You can run the file test_final_from_paper.m to run exactly the same experiments as in our paper. 


It also contains the following codes for the same problem: 

- from the paper 'F.-X. Orbandexivry, Yu. Nesterov, and P. Van Dooren, Nearest stable system using successive convex approximations, Automatica, 49 (2013), pp. 1195-1203' that was kindly made available to us by F.-X. Orbandexivry. Of course this code remains the property of these authors and if using this code, please cite the paper above. 

- from Michael Overton based on his toolbox HANSO; see http://www.cs.nyu.edu/overton/software/hanso/. Similarly as above, the code remains the property of M. Overton and co-authors should be cited properly if their code is used. 

- from Tim Mitchell, a method based on GRANSO that solves non-smooth optimization problems; see the paper http://dx.doi.org/10.1080/10556788.2016.1208749. Note that comparisons with this algorithm are not provided in our paper since Tim provided us with this algorithm later on. For large matrices, our code still works best (for small matrices, it depends to which local minima each algorithm converges to; cf. the discussion in the paper). 
To run this code, you need to download the GRANSO toolbox http://timmitchell.com/software/GRANSO/ 



