! ALF_x module for computing fully normalized associated Legendre functions
! using extended range arithmatic method. this software provided by Fukushima (2012)
! small change applied to adopted with our software     
!    Fukushima, T. Numerical computation of spherical harmonics of arbitrary degree 
!    and order by extending exponent of floating point numbers. J Geod 86, 271–285 (2012) 
    module ALF_x
    
    use constant

    implicit none

    contains
  
    subroutine alfx (m, nx, t, xi, eta , ps, ips, p)

    implicit none
    integer, intent (in) :: m, nx, ips
    real(kind=Dp), intent(in):: ps, t
    real(kind=dp) p1n2, p1, anm, bnm, p1n1, temp, pin2, pin1, pin, p1n
    integer::  n, ip1n2, ip1n1, n1, ip1n
    real(kind=dp),  dimension (:),  allocatable , intent(in) :: xi, eta
    real(kind=dp),  dimension (:),  allocatable, intent(inout) ::  p
 
    p1n2=ps;ip1n2=ips;p1=x2f(p1n2,ip1n2);p(m)=p1

    if(m+1.gt.nx) return

    anm=xi(2*m+3);

    p1n1=anm*t*ps;

    ip1n1=ips;

    call xnorm(p1n1,ip1n1)

    p1=x2f(p1n1,ip1n1);

    p(m+1)=p1

    do n=m+2,nx
    
        temp=xi(2*n+1)*eta(n+m)*eta(n-m);
        anm=temp*xi(2*n-1)
        bnm=temp*eta(2*n-3)*xi(n+m-1)*xi(n-m-1)
        call xlsum2(anm*t, p1n1 , ip1n1, -bnm, p1n2, ip1n2, p1n, ip1n);

         p1=x2f(p1n,ip1n)
         p(n)=p1

        pin2=pin1;
        pin1=pin;
        p1n2=p1n1;
        p1n1=p1n;

        ip1n2=ip1n1;
        ip1n1=ip1n;

        if(ip1n2.eq.0) goto 1

    enddo

    return
1   continue
   
    n1=n+1
    
     do n=n1,nx
        temp=xi(2*n+1)*eta(n+m)*eta(n-m);
        anm=temp*xi(2*n-1)
        bnm=temp*eta(2*n-3)*xi(n+m-1)*xi(n-m-1)
         p(n) = anm*t*p(n - 1) - bnm *p(n - 2);
    
    enddo

    end subroutine alfx

    !==================================================================
    subroutine alfsx(mx, d, u1, ps,ips)
    implicit none

    integer mx,ix1,m
    real(kind=dp), allocatable::  d(:) , ps(:)
    integer, allocatable:: ips(:)
    real(kind=dp)::  u1
    real(kind=dp) x1, dm

    ps(0)=1.0_dp;
    ips(0)=0

    x1=sqrt(3.0_dp)*u1;
    ix1=0;
    ps(1)=x1;
    ips(1)=ix1

    dm=d(2)
    x1=(dm*u1)*x1;
    call xnorm(x1,ix1);
    ps(2)=x1;
    ips(2)=ix1

    do m=3,mx
        dm=d(m)
        x1=(dm*u1)*x1;
        call xnorm(x1,ix1);
        ps(m)=x1;
        ips(m)=ix1
    enddo
    return;
    end subroutine alfsx
    !==========================================================================

    subroutine xnorm(x,ix)
    integer, intent(inout):: ix
    real(kind=dp), intent(inout):: x
    real(kind=dp):: y
    y= abs(X)
    if(x.eq.0.0_dp) then
        ix=0
    elseif(y.LT.BIGSI) then 
        x=x*BIG;
        ix=ix-1
    elseif(y.GE.BIGS) then 
        x=x*BIGI;
        ix=ix+1
    endif
    return
    end subroutine xnorm
    !====================================================================
   subroutine xlsum2(a,x,ix,b,y,iy,z,iz)

    real(kind=dp), intent(in) :: a, x, b, y
    real(kind=dp):: q
    real(kind=dp), intent(out):: z
    integer, intent(in):: ix, iy
    integer, intent(out):: iz
    integer iq, id

    z=a*x; iz=ix
    q=b*y; iq=iy
    id=iz-iq
    if(id.eq.0) then
        z=z+q
    elseif(id.eq.1) then
        z=z+q*BIGI
    elseif(id.eq.-1) then
        z=q+z*BIGI; iz=iq
    elseif(id.lt.-1) then
        z=q; iz=iq
    endif
    call xnorm(z,iz)
    return
    end  subroutine xlsum2
    !==========================================================================

    real(kind=dp) function x2f(x,ix)
    real(kind=dp), intent(in):: x
    integer, intent(in) :: ix

    if(ix.eq.0) then
        x2f=x
    elseif(ix.eq.-1) then
        x2f=x*BIGI
    elseif(ix.eq.1) then
        x2f=x*BIG
    elseif(ix.lt.0) then
        x2f=0.0_dp
    elseif(x.lt.0) then
        x2f=-BIG
    else
        x2f=BIG
    endif
    return
    end function x2f
  
    
    end module ALF_x

