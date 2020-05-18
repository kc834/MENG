# MENG

Kevin Chaudhari
Title: Comparison of Bi-Lanczos and Arnoldi Methods
M.Eng. Project, in partial fulfilment of the Cornell CS M.Eng. Degree. 

arnoldi.m, lanczos.m, bilanczos.m are stand-alone versions of their respective algorithms and are located in the standalone folder.

expv.m is an Arnoldi solver for exp(-At)v developed as part of the EXPOKIT library.

expvB.m is a Bi-Lanczos solver for u^{\dagger}\exp(-At)v and is inspired by the aforementioned Arnoldi solver.

compare.m runs both expv and expvB on a given u, A, v, t, tolerance, and dimension. The output is the relative error between the two implementations.

utfAv.m runs expv and expvB on a user-provided A_matrix.dat, u.dat, and v.dat file. This is the main function. 

test1_real.m runs expv and expvB on a randomly generated 800-by-800 sparse real matrix and randomly generated 800-by-1 real starting vectors.

test2_complex.m runs expv and expvB on a randomly generated 800-by-800 sparse complex matrix and randomly generated 800-by-1 complex starting vectors.

Miscellaneous: The workspace folder contains .mat files with some example values of A, u, v and preset parameters for additional testing. importA.m is used to extract A_matrix.dat but can be replaced by the user depending on how the matrix is stored.

References:
1. Derived from lecture notes of David Bindel, CS 6210, Fall 2019, Cornell
2. High-level algorithm from PhD Thesis of Axel Facius, July 2000, Universitat Karlsruhe
3. Rational approximation derived with insight from "Efficient Solution of Parabolic Equations by Krylov Approximation Methods" by E. Gallopoulos and Y. Saad
4. Insight for preprocessing from "Analysis of Some Krylov Subspace Approximations to the Matrix Exponential Operator" by Y. Saad
5. Time-stepping protocol from "Expokit: A Software Package for Computing Matrix Exponentials" by Roger B. Sidje