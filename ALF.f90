!module ALF
!    module for computation fully normalized associated Legerndre function
!    programmer : Mehdi Goli
!    Goli@shahroodut.ac.ir
    
    module ALF

    use constant

    implicit none

    contains

    !-------------------------------------------------------

    subroutine ini_ALF (nmax, xi, eta, d, logd)
    
    ! precompute some arrays to fast evaluation an, bn terms
    ! based on Fukushima (2012) code
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    ! input: nmax : maximum degree of expantion
    ! output: xi, eta, d, logd
    
    implicit none
    integer::nmax
    integer:: m
    real(kind=dp),  dimension (:),  allocatable:: xi, eta, d
    real(8),  dimension (:),  allocatable, optional:: logd
    
    allocate(eta(3*nmax+3), xi(-1:3*nmax+3),d(nmax))
    
    xi(-1:0)=0d0;

    do m=1,3*nmax+3
        xi(m)=sqrt(dble(m))
        eta(m)=1.0_dp/xi(m)
    enddo

    do m=1,nmax
        d(m)=xi(2*m+1)*eta(2*m)
    enddo
    if (present(logd)) then
        allocate(logd(nmax))
        logd=log(d)
    endif
    
    end subroutine ini_ALF

    !-------------------------------------------------------
   
    subroutine alfs(mx, d, u, ps)
    implicit none

    !  sectorial fully normalized ALF pmm(x)
    ! if pmm<2d-960 the pmm is replaced by zero    
    ! input: mx : maximum degree of expantion
    !        d    : precomputed values of sqrt(2n / 2n+1) 
    !        u    : sqrt(1-x^2)
    !        
    ! output: ps(0:mx) :sectorial fully normalized ALF pmm(x)
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    
    integer mx,m
    real(kind=dp), allocatable::  d(:) , ps(:)

    real(kind=dp)::  u
    real(kind=dp) x1, dm

    ps(0)=1d0;

    x1=sqrt(3d0)*u;

    ps(1)=x1;

    dm=d(2)
    x1=(dm*u)*x1;
    ps(2)=x1;

    do m=3,mx
        dm=d(m)
        x1=(dm*u)*x1;
        ps(m)=x1;
        
        if (x1<bigi) then
            ps(m:mx)=0d0
            return
        endif
        
    enddo
    return;
    end

    !-------------------------------------------------------

    subroutine alfs_log(mx, d, u, logd, ps, log_ps)
  

    !   lograithm of sectorial fully normalized ALF pmm(x)
    !   
    ! input: mx : maximum degree of expantion
    !        d    : precomputed values of sqrt(2n / 2n+1) 
    !        u    : sqrt(1-x^2)
    !       logd  : precomputed values of log(d)
    ! output: ps(0:mx) :sectorial fully normalized ALF pmm(x)
    !       log_ps: log of sectorial fully normalized ALF pmm(x)
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    
      implicit none
    integer mx,m
    real(kind=dp), allocatable::  d(:) , ps(:), log_ps(:)
    
    real(8), allocatable::  logd(:)

    real(kind=dp):: u
    
    real(8):: log_u
    
    real(kind=dp):: x1, dm, up 
 
    ps(0)=1d0;
   
    log_ps(0)=0d0

    x1=sqrt(3d0)*u;

    ps(1)=x1;
    
    up=u
    
    m=0
    do while (up>0.01)
        up=up/2d0
        m=m+1
    enddo
    
    log_u=log(up)+m*log(2d0)

    log_ps(1)=logsqrt3+log_u
    
    do m=2,mx
        
        dm=d(m)
        
        x1=(dm*u)*x1;
        ps(m)=x1;
       
        log_ps(m)=log_ps(m-1)+ logd(m)

    enddo
    
    do m=2,mx
        log_ps(m)=log_ps(m)+(m-1)*log_u
    enddo
    
    return;
    end

    !-------------------------------------------------------

    subroutine alfs_ratio(mx, d, u, ps, ips)
    implicit none

    !  sectorial fully normalized ALF pmm(x) using ERA
    !   
    ! input: mx : maximum degree of expantion
    !        d    : precomputed values of sqrt(2n / 2n+1) 
    !        u    : sqrt(1-x^2)
    !  
    ! output: ps(0:mx) : priciple part of  pmm(x)
    !         ips(0:mx): auxiliary index of pmm(x)
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    
    integer mx,m, ix
    real(kind=dp), allocatable::  d(:) , ps(:)
    integer, allocatable::  ips(:)
    real(kind=dp)::  u
    real(kind=dp):: x1, dm

    ix=0;

    ps(0)=1d0;
    ips(0)=0

    x1=sqrt(3d0)*u;

    ps(1)=x1;
    ips(1)=0

    dm=d(2)
    x1=(dm*u)*x1;
    ps(2)=x1;
    ips(2)=0

    do m=3,mx

        dm=d(m)
        x1=(dm*u)*x1;

        call tiny_f(x1, ix)

        ps(m)=x1;
        ips(m)=ix

    enddo

    return;
    end

    
    !-------------------------------------------------------

    subroutine alf_mid(m, nmax, t, u, xi, eta , ps, pm)

    !  fully normalized ALF, P_(m:nx,m)(x) using midway method
    !   
    ! input: m  : order m<=mx
    !        mx : maximum degree of expantion
    !        t  : x = t = cos(theta), theta: co-latitude
    !        u  : sqrt(1-t^2)
    !        xi : precomputed values of sqrt(n)
    !        eta: precomputed values of 1/sqrt(n)
    !        ps : P_mm sectorial ALF
    !        pm : fully normalized ALF previous order, i.e., P_(m-1:nx,m-1)
    !  
    ! output: pm(0:mx) : fully normalized ALF  
    !         note that pm(0:m-1)=0
    !         pm is input/output argument.
    !      
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    implicit none

    real(kind=dp) t, u,  ps, temp, anm, bnm, pm_n, pm_np1
    integer:: m, nmax, nmin, n
    real(kind=dp),  dimension (:),  allocatable:: xi, eta
    real(kind=dp),  dimension (:),  allocatable:: pm

    if (m==nmax) then

        pm(nmax)=ps
        return

    endif

    if (ps >= bigi ) then

        pm(m)=ps

        anm=xi(2*m+3);

        pm(m+1)=anm*t*ps;

        nmin=m

    else

        do n=m, nmax
            if ( abs(pm(n)) > bigi) exit
        enddo

        nmin=n

        if (n>nmax-3) return

        anm=xi(n-m+1)*eta(n+m);
        bnm=   xi(2*n+1)*xi(n+m-1) *eta(2*n-1)*eta(n+m);
        pm_n= -anm*t/u*pm(n)+ 1d0/u*bnm*pm(n-1);

        n=n+1

        anm=xi(n-m+1)*eta(n+m);
        bnm=   xi(2*n+1)*xi(n+m-1) *eta(2*n-1)*eta(n+m);
        pm_np1= -anm*t/u*pm(n)+ 1d0/u*bnm*pm(n-1);

        pm(n-1)=pm_n
        pm(n)  =pm_np1

    endif

    do n=nmin+2, nmax

        temp=xi(2*n+1)*eta(n+m)*eta(n-m);
        anm=xi(2*n-1)*temp
        bnm=eta(2*n-3)*xi(n+m-1)*xi(n-m-1)*temp
        pm(n)=anm*t*pm(n-1)-bnm*pm(n-2);

    enddo

    end

    !-------------------------------------------------------

    subroutine alf_ratio_log( m, L, t, xi , eta, ps , log_ps, pn )

    !  fully normalized ALF, P_(m:nx,m)(x) using successive ratio method using logarithm
    !   
    ! input: m  : order m<=mx
    !        mx : maximum degree of expantion
    !        t  : x = t = cos(theta), theta: co-latitude
    !        u  : sqrt(1-t^2)
    !        xi : precomputed values of sqrt(n)
    !        eta: precomputed values of 1/sqrt(n)
    !        ps : P_mm sectorial ALF
    !        log_ps : log of P_mm
    !  
    ! output: pn(0:mx) : fully normalized ALF  
    !      
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    
    integer, intent (in):: m, L

    real(8), dimension(:), allocatable, intent(in):: xi, eta!, pch

    real(8), dimension(:), allocatable, intent(inout):: pn

    real(8), intent(in):: log_ps , t , PS

    real(8):: temp,  anm, bnm , d1, d2, d

    integer:: n , n0 , j, nlog

    real(8):: log0
  
    if (log_ps > log_tol) then

        pn(m)=ps

        if (m==L) return

        anm=xi(2*m+3);

        pn(m+1)=anm*t*ps;

        do n=m+2,L

            temp=xi(2*n+1)*eta(n+m)*eta(n-m);
            anm=temp*xi(2*n-1)
            bnm=temp*eta(2*n-3)*xi(n+m-1)*xi(n-m-1)
            pn(n)=anm*t*pn(n-1)-bnm*pn(n-2);

        enddo

    else

        log0=log_ps;

        n=m+1;

        anm=xi(2*n-1)*xi(2*n+1)*eta(n-m)*eta(n+m);

        d2 = 1d0

        d1 =  anm*t

        do n=m+2, L

                temp=xi(2*n+1)*eta(n+m)*eta(n-m);

                anm =xi(2*n-1)*temp

                bnm= eta(2*n-3)*xi(n+m-1)*xi(n-m-1) *temp

                d= anm*t*d1 - bnm*d2 ;

                d2=d1

                d1=d

                if (d>bg) then
                    d=d*bgi
                    d1=d1*bgi
                    d2=d2*bgi
                    log0=log0+log_bg
                endif
                
                if (log0>log_tol) exit
              
        enddo
        
        if (n>L-1) return

        pn(m:n-2)=0d0;

        log0=log0+log(d)

        pn(n)= exp(log0) ;

        pn(n-1)=  pn(n)*d2/d1

        n0=n+1  ;

        do n=n0,L

            temp=xi(2*n+1)*eta(n+m)*eta(n-m);
            anm=temp*xi(2*n-1)
            bnm=temp*eta(2*n-3)*xi(n+m-1)*xi(n-m-1)
            pn(n)=anm*t*pn(n-1)-bnm*pn(n-2);

        enddo
       
    endif

    end


    !-------------------------------------------------------

    subroutine alf_ratio_x( m, L, t, xi , eta, ps , ips, pn)

    
    !  fully normalized ALF, P_(m:nx,m)(x) using succesive ratio method using ERA
    !   
    ! input: m  : order m<=mx
    !        mx : maximum degree of expantion
    !        t  : x = t = cos(theta), theta: co-latitude
    !        u  : sqrt(1-t^2)
    !        xi : precomputed values of sqrt(n)
    !        eta: precomputed values of 1/sqrt(n)
    !        ps/ips : principal and auxiliary part of P_mm
    !        
    ! output: pm(0:mx) : fully normalized ALF  
    !         note that pm(0:m-1)=0
    !      
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    implicit none

    integer, intent (in)::  m, L, ips

    real(8), allocatable, intent(in):: xi(:), eta(:)
    real(8), allocatable, intent(inout):: pn(:)

    real(8), intent(in):: t , PS

    real(8)::  temp,  anm, bnm

    integer:: n, ix , n0, j

    real(8):: d1, d2, d

    if (ips == 0) then

        pn(m)=ps

        if (ips==1) pn(m) = pn(m) * BiGsI

        if (m==L) return

        anm=xi(2*m+3);

        pn(m+1)=anm*t*pn(m);

        do n=m+2,L

            temp=xi(2*n+1)*eta(n+m)*eta(n-m);
            anm=temp*xi(2*n-1)
            bnm=temp*eta(2*n-3)*xi(n+m-1)*xi(n-m-1)

            pn(n) = anm*t*pn(n - 1) - bnm * pn(n - 2);

        enddo

    else

        if (m==L) return

        n=m+1;

        anm=xi(2*n-1)*xi(2*n+1)*eta(n-m)*eta(n+m);

        d2 = ps

        d1= ps *anm*t

        ix=ips

        do n=m+2, L

            temp=xi(2*n+1)*eta(n+m)*eta(n-m);

            anm =xi(2*n-1)*temp

            bnm= eta(2*n-3)*xi(n+m-1)*xi(n-m-1) *temp

            d= anm*t*d1 - bnm*d2 ;

            d2=d1

            d1=d

            if (ix == 1)then

                pn(n-1)=d2* bigi
                pn(n)=d*bigi
                if (pn(n-1) > bigi) exit

            endif

            if (d>=BiGs) then

                d2= d2*BiGI
                d1= d1*BiGI
                d= d*BiGI
                ix=ix-1

            endif

        enddo

        if (n>L)return

        n0=n+1;

        do n=n0,L

            temp=xi(2*n+1)*eta(n+m)*eta(n-m);
            anm=temp*xi(2*n-1)
            bnm=temp*eta(2*n-3)*xi(n+m-1)*xi(n-m-1)

            pn(n) = anm*t*pn(n - 1) - bnm *pn(n - 2);

        enddo

    endif

    end
   
    !-------------------------------------------------------

    pure subroutine tiny_f(x, ix)
    real(kind=Dp), intent(inout):: x
    integer, intent(inout):: ix
    
    ! rescale the positive real number x when x<BGsI 
    ! input/output :
    !  x, ix:  principal and auxiliary part of P_mm
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    if (x<=BiGsI) then
        x=x*BiG
        ix=ix+1
    endif

    end

    !-------------------------------------------------------
    pure subroutine big_f(x, ix)
    real(kind=Dp), intent(inout):: x
    integer, intent(inout):: ix

    
    ! rescale the positive real number x when x>BGs 
    ! input/output :
    !  x, ix:  principal and auxiliary part of P_mm
    !    programmer : Mehdi Goli
    !    Goli@shahroodut.ac.ir
    if (x>=BiGs) then
        x=x*BiGI
        ix=ix-1
    endif

    end
    
    !---------------------------------------------------------------------

    end module ALF
