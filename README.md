# DASKR 

## Preamble

This repository is a mirror of the netlib code of [daskr](https://www.netlib.org/ode/). It is stored here in 
the hope to be modernized some day. 
At the moment, as compared to the original code, it has been restructured to be compatible with **fpm** build system.

TODO list
- [x] Use fpm as build system
- [ ] Remove implicite typing
- [ ] Remove goto's
- [ ] Create explicit interfaces to pass residual, jacobian and krilov solver.
- [ ] Split direct and krylov solvers.
- [ ] Replace linpack linear solver for Lapack
- [ ] Create a modern documentation site

Since the code is a the the moment F77 fixed-form it should be compiled with the option `-std=legacy`. 

```cmd
fpm build --flag "-std=legacy"
```

## Introduction

       DASKR Package: DAE Solver with Krylov Methods and Rootfinding
                      Version of 8 June 2011

            P. N. Brown, A. C. Hindmarsh, and L. R. Petzold 


DASKR is a solver for systems of differential-algebraic equations (DAEs).  
It includes options for both direct and iterative (Krylov) methods for the
solution of the linear systems arising at each (implicit) time step.
DASKR is a variant of the DASPK package [1].  In addition to all the
capabilities of DASPK, DASKR includes the ability to find the roots of a
given set of functions while integrating the DAE system.

In contrast to the older DASSL package, DASKR includes a procedure
for calculating consistent initial conditions for a large class of
problems (which includes semi-explicit index-1 systems) [2].  This
procedure includes options for inequality constraints on selected
components.  The package also includes an option to omit the algebraic
components from the local error control.

Along with the solver itself, the DASKR package includes five example
programs and a set of general-purpose preconditioner files.  These are
described in more detail below.  The package includes separate single
and double precision versions of all source files.


Package Contents
----------------

1. The DASKR package is being distributed in the form of a tar file,
which expands to a directory, DASKR.  The DASKR directory contains
this README file and three subdirectories.

2. The subdirectory DASKR/solver contains the following source files:

     ddaskr.f = main solver source, double precision
     dlinpk.f = required LINPACK/BLAS routines, double precision
     daux.f   = machine constant and error handler, double precision

     sdaskr.f = main solver source, single precision
     slinpk.f = required LINPACK/BLAS routines, single precision
     saux.f   = machine constant and error handler, single precision

The LINPACK/BLAS files are provided for the sake of completeness, but
the target machine/system of interest may already have optimized versions
of these routines, which should be used instead for the sake of efficiency.

A complete usage document is included as the initial prologue of each
source file ddaskr.f/sdaskr.f.

3. The subdirectory DASKR/examples contains the following files:

     dkrdem.f   = small demonstration program, double precision, which
                  solves two small problems with self-checking logic
     dheat.f    = heat equation example program, double precision, using
                  the Krylov option with banded preconditioner
     dheatilu.f = heat equation example program, double precision, using
                  the Krylov option with sparse ILU preconditioner
     dweb.f     = food web system example program, double precision, using
                  the direct (band) option, and the Krylov option with a
                  product preconditioner based on reaction-transport form
     dwebilu.f  = food web system example program, double precision, using
                  the Krylov option with sparse ILU preconditioner
     makeddem   = Make-file to compile and load dkrdem program
     makedh     = Make-file to compile and load dheat program
     makedhilu  = Make-file to compile and load dheatilu program
     makedw     = Make-file to compile and load dweb program
     makedwilu  = Make-file to compile and load dwebilu program

     skrdem.f   = small demonstration program, single precision, which
                  solves two small problems with self-checking logic
     sheat.f    = heat equation example program, single precision, using
                  the Krylov option with banded preconditioner
     sheatilu.f = heat equation example program, single precision, using
                  the Krylov option with sparse ILU preconditioner
     sweb.f     = food web system example program, single precision, using
                  the direct (band) option, and the Krylov option with a
                  product preconditioner based on reaction-transport form
     swebilu.f  = food web system example program, single precision, using
                  the Krylov option with sparse ILU preconditioner
     makesdem   = Make-file to compile and load skrdem program
     makesh     = Make-file to compile and load sheat program
     makeshilu  = Make-file to compile and load sheatilu program
     makesw     = Make-file to compile and load sweb program
     makeswilu  = Make-file to compile and load swebilu program

     dkrdem.out   = output from dheat.f program
     dheat.out    = output from dheat.f program
     dheatilu.out = output from dheatilu.f program
     dweb.out     = output from dweb.f program
     dwebilu.out  = output from dwebilu.f program

Except for dkrdem/skrdem, all of these examples make use of preconditioner
routines, provided in the files described below, when calling DASKR with
the Krylov method option.  All are heavily commented and intended for use
as models for real user applications of DASKR.

4. The subdirectory DASKR/preconds contains the following source files:

     dbanpre.f = preconditioner for banded problems, double precision
     drbdpre.f = routines for a reaction-based block-diagonal precond-
                 itioner for DAEs arising from a reaction-transport
                 system, without block-grouping, double precision
     drbgpre.f = routines for a reaction-based block-diagonal precond-
                 itioner for DAEs arising from a reaction-transport
                 system, with block-grouping, double precision
     dilupre.f = preconditioner routines for sparse ILU preconditioning,
                 double precision, intended for general DAE problems.
     dsparsk.f = routines from SPARSKIT, Y. Saad (Univ. of Minnesota),
                 double precision, used by dilupre.f

     sbanpre.f = preconditioner for banded problems, single precision
     srbdpre.f = routines for a reaction-based block-diagonal precond-
                 itioner for DAEs arising from a reaction-transport
                 system, without block-grouping, single precision
     srbgpre.f = routines for a reaction-based block-diagonal precond-
                 itioner for DAEs arising from a reaction-transport
                 system, with block-grouping, single precision
     silupre.f = preconditioner routines for sparse ILU preconditioning,
                 single precision, intended for general DAE problems.
     ssparsk.f = routines from SPARSKIT, Y. Saad (Univ. of Minnesota),
                 single precision, used by silupre.f

At the beginning of each of these source files is a prologue which documents
the contents of the file and the usage of the routines in it.


Installation and Usage Notes
----------------------------

1. The single and double precision versions of the SPARSKIT subset, in
the files ssparsk.f and dsparsk.f in DASKR/preconds, cannot be
installed together as a single library, because of name duplications
among the individual routines.  If such a combined installation is
desired, first change the names of all precision-dependent subroutines
in either ssparsk.f or dsparsk.f, so as to have unique names across
the two files.  Then make the same name changes in the corresponding
file silupre.f or dilupre.f which calls those SPARSKIT routines.

2. The five example problems can be compiled and loaded using the
appropriate make-file in DASKR/examples.  First check the compiler and
flags in the make-file, however.  Each make-file finds the required
solver and preconditioner source or object files in ../solver and
../preconds.  The output files for all five examples, as run in double
precision on a Sun Sparc-10 Workstation, are provided in the files
d*.out.  The results on other systems may differ slightly.

3. The example programs that use ILU preconditioning, if altered to
use the option JACOUT = 1, write an additional output file containing
the initial Jacobian and residual vector in Boeing-Harwell format.  A
map of the Jacobian as a postscript file can then be generated using
the SPARSKIT routines readmt and pspltm.

4. Users of DASKR with a sparse ILU preconditioner are encouraged to
experiment with the many options available in the SPARSKIT package,
as accessed through the dilupre.f/silupre.f modules provided.  See the
prologue in those files, and the two examples that use ILU.

5. Important Note: The use of single precision on a 32-bit machine is
strongly discouraged.  The amplification of roundoff errors by the
integration and linear system solution algorithms is such that double
precision is generally required on short-wordlength machines.  In
particular, the small demonstration program and the food web example
programs provided here do not run on a 32-bit machine in single
precision, with even the moderate tolerance values of 1.0e-5.


References
----------
[1] P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov
    Methods in the Solution of Large-Scale Differential-Algebraic
    Systems, SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488.

[2] P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent
    Initial Condition Calculation for Differential-Algebraic
    Systems, SIAM J. Sci. Comp. 19 (1998), pp. 1495-1512.

---------------------------------------------------------------------------

Work performed under the auspices of the U.S. Department of Energy by
Lawrence Livermore National Laboratory under contract number W-7405-Eng-48.

Copyright (c) 2002, The Regents of the University of California. 
UCRL-CODE-2002-058
All rights reserved. 
This file is part of DASKR.
For details, see DASKR/LICENSE
