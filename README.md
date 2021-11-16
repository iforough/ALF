# ALF
Computation of the Associated Legendre Functions up to ultra high degree/order

Version 1.0 - 2021

Mehdi Goli - Goli@shahroodut.ac.ir
Ismael Foroughi - foroughi.ismael@gmail.com

ALF_UHD computes the fully normalized associated Legendre function (fnALF)
using fixed-order increasing degree recursive formula for arbitrary 
degree n, and order m (Holmes and Featherstone, 2002):

P_(nm) (x)= anm * x * P_(n-1,m) -
            bnm * P_(n-2,m) 
n,m are degree and order
x = cos(theta), theta is co-latitude
note that:  theta /= 0, \pi

The code should compile on any Linux/WINDOWS system using Intel or gfortran 
compilers. We suggest the gfortran compiler under optimization command
'-O3 -march=native' to get maximum performance. 

ALF_UHD involved three modules
    - constants
    module to define parameters 

    - ALf
    module to compute the fnALF using successive ratio of fnALF.

    This module contains:
    
    -- subroutine ini_ALF
           ini_ALF (pre)computes some arrays to fast evaluation 
           an, bn terms 
     
   -- subroutine alfs 
           alfs computes sectorial fnALF pmm(x) > 2^-960 (tiny values 
           of pmm(x) is replaced by zero).
     
   -- subroutine alfs_log 
           alfs computes sectorial fnALF pmm(x) > 2^-960 and log(Pmm).
     
   -- subroutine alfs_ratio
           alfs_ratio computes the principle and auxiliary parts of 
           sectorial fnALF pmm(x) using extended range arithmatic (ERA).
     
   -- subroutien midway
           midway compute efficiently fnALF using midway method. The seed 
           values were evaluated using previous order when underflow error
           accures.
    
    -- subroutine alf_ratio_log
           alf_ratio_log computes fnALF using successive ratio of fnALF. 
           The underflow problem handled using logarithm.
    
    -- subroutine alf_ratio_x
           alf_ratio_x computes fnALF using successive ratio of fnALF. 
           The underflow problem handled using ERA.
      
    -- subroutine tiny_f(x, ix)
           tiny_f changes principle and auxiliary parts of positive real
           number x when x < BGsI. BGSI=2^960 is a big number.  
        
    -- subroutine big_f(x, ix)
           big_f changes principle and auxiliary parts of positive real
           number x when x > BGs. BGS=2^-960 is a tiny number. 

    - ALf_x 
    module to compute the fnALF using EAR method.
    ALF_x was developed by Fukushima (2012). We changed it a bit to adapt
    to our software.

examples:

    - time_alf  program compares the run-time of fnALf computation for
                 different Nmax. The following method applied in time_alf:
                          1- successive ratio with logarithm
                          2- successive ratio with ERA
                          3- midway method
                          4- ERA using Fukushima's code
                   
    - time_alf_optimize program gives the run-time of optimized fnALf 
                 computation for different Nmax. The following method 
                 applied in time_alf_optimize:
                          1- successive ratio with logarithm
                          2- successive ratio with ERA
                          3- midway method
                          4- ERA using Fukushima's code
                   
