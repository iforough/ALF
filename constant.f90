module constant
    
    integer, parameter :: sp = selected_real_kind(6, 33)
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(kind=dp), parameter:: pi=3.141592653589793d0!23846264338327950288_dp
    real(kind=dp), parameter:: dr=pi/180.0_dp !degree2radius
    integer,parameter:: IND=960
    real*8, parameter:: BIG=2.d0**(IND )
    real*8, parameter:: BIGI=2.d0**(-IND )
    real*8, parameter:: BIGS=2.d0**(IND/2)
    real*8, parameter:: BIGSI=2.d0**(-IND/2)
    real*8, parameter:: log_tol = -665.421293337547497040542836599850
    real*8, parameter:: log_big = log(big)
    
    real*8, parameter:: BG=1d30
    real*8, parameter:: BGI=1d-30
    real*8, parameter:: log_bg=log(BG)
     
    real(8):: logsqrt3= 0.54930614433405484569d0
    real(8):: log2 =0.693147180559945309417232121458177
      
    end module constant 