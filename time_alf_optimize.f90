    !gfortran Alf.f90 ALfx.f90 constant.f90 time_alf_optimize.f90 -O3 -march=native

    use alf
    use alf_x
    use constant

    implicit none

    integer:: L, m, n, i, k, j, n0
    real(8), allocatable:: ps(:), pm(:), xi(:), eta(:), d(:), log_ps(:), logd(:)
    integer, allocatable:: ips(:)
    real(8):: theta, u, t, t1, t2, t3, t4, step


    open(1452, file='time_analysis_opt.txt')

    print*, 'cpu time for optimize computation fnALF for uniform grid over theta in [0,89] : '
    print*
    print*, 'Nmax  lat_number   ratio_ERA   ratio_log  Fukushima midway'

    step=1d-6

    do i= 4, 16

        L=2**i

        step=step*4

        if (step>1) step=1d0
        k=89d0/step
 
        write(1452, '(i6,2x, i10, 2x    )',advance='no') L, k

        write(*, '(i6,2x, i10 , 2x)',advance='no') L, k

        !======================================================= ERA_ratio

        if (allocated(ps)) deallocate(ps)
        if (allocated(pm)) deallocate(pm)
        if (allocated(ips)) deallocate(ips)
        if (allocated(d)) deallocate(d)
        allocate( ps(0:L), ips(0:L), pm(0:L))
        if (allocated(xi)) deallocate(xi)
        if (allocated(eta)) deallocate(eta)
        if (allocated(d)) deallocate(d)

        call ini_ALF (L, xi, eta, d)

        call cpu_time(t1)

        theta=0d0

        do j=1, k

            theta=theta+step

            u=sin(dr*theta)
            t=cos(dr*theta)

            call alfs_ratio(L, d, u, ps , ips)

            n0=min(int(L*u+200),L)

            do m=0, n0!L

                call alf_ratio_x( m, L, t,  xi, eta , ps(m), ips(m), pm)

            end do

        enddo

        call cpu_time(t2)

        write(1452, '(  2x, f8.3 )',advance='no')   t2-t1

        write(*, '( 2x, f8.3  )',advance='no')  t2-t1


        !==================================================log_ratio

        if (allocated(ps)) deallocate(ps)
        if (allocated(pm)) deallocate(pm)
        if (allocated(log_ps)) deallocate(log_ps)
        if (allocated(d)) deallocate(d)
        allocate( ps(0:L), log_ps(0:L), pm(0:L))
        if (allocated(xi)) deallocate(xi)
        if (allocated(eta)) deallocate(eta)
        if (allocated(d)) deallocate(d)
        if (allocated(logd)) deallocate(logd)

        call ini_ALF (L, xi, eta, d, logd=logd)

        call cpu_time(t1)

        theta=0d0

        do j=1, k

            theta=theta+step

            u=sin(dr*theta)
            t=cos(dr*theta)

            call  alfs_log(L, d, u, logd, ps , log_ps)

            n0=min(int(L*u+200),L)

            do m=0, n0!L

                call alf_ratio_log( m, L, t,  xi, eta , ps(m), log_ps(m), pm)

            end do

        enddo

        call cpu_time(t2)

        write(1452, '(  2x, f8.3  )',advance='no')   t2-t1

        write(*, '( 2x, f8.3 )',advance='no')  t2-t1


        !======================================= Fukushima
        if (allocated(ps)) deallocate(ps)
        if (allocated(pm)) deallocate(pm)
        if (allocated(d)) deallocate(d)
        if (allocated(ips)) deallocate(ips)

        allocate( ps(0:L), ips(0:L), pm(0:L))
        if (allocated(xi)) deallocate(xi)
        if (allocated(eta)) deallocate(eta)
        if (allocated(d)) deallocate(d)

        call ini_ALF (L, xi, eta, d)

        call cpu_time(t1)

        theta=0d0

        do j=1, k

            theta=theta+step

            u=sin(dr*theta)
            t=cos(dr*theta)

            call alfsx(L, d, u, ps, ips)

            n0=min(int(L*u+200),L)

            do m=0, n0!L

                call alfx(m, L, t,  xi, eta , ps(m), ips(m), pm)

            end do

        enddo

        call cpu_time(t2)

        write(1452, '(  f8.3 )',advance='no')    t2-t1

        write(*, '(  f8.3 )',advance='no')     t2-t1

        !=======================================================midway

        if (allocated(ps)) deallocate(ps)
        if (allocated(pm)) deallocate(pm)
        if (allocated(log_ps)) deallocate(log_ps)
        allocate( ps(0:L), log_ps(0:L), pm(0:L))
        if (allocated(xi)) deallocate(xi)
        if (allocated(eta)) deallocate(eta)
        if (allocated(d)) deallocate(d)
        call ini_ALF (L, xi, eta, d)

        call cpu_time(t1)

        theta=0d0
        
        do j=1, k

            theta=theta+step

            u=sin(dr*theta)
            t=cos(dr*theta)

            call alfs(L, d, u, ps )

            n0=min(int(L*u+200),L)

            do m=0, n0!L

                call alf_mid( m, L, t, u, xi, eta , ps(m) , pm)

            end do

        enddo

        call cpu_time(t2)

        write(1452, '(  2x, f8.3 )')   t2-t1

        write(*, '( 2x, f8.3 )')  t2-t1
 
    enddo

    end

