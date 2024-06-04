module convection_module
contains

    subroutine initvalue(un, u, ua, c, c0, c1, x, k, nx, signal_type)
        integer :: signal_type, nx, i
        real*16 :: k
        real*16, dimension(nx) :: x 
        real*16, dimension(0:nx+1) :: un, u, ua
        real*16 :: c0, c1, c

        if (signal_type == 1) then
            un = c0 + c1*sin(k*x)
            u = un
        elseif (signal_type == 2) then
            do i = 1, nx
                if ((x(i)>=0.4) .and. (x(i)<=0.7)) then
                    un(i) = 20.0*x(i)/3.0-5.0/3.0
                else
                    un(i) = 1.0
                end if
                u = un
            end do

            do i = 1, nx
                if ((x(i)>=0.4+c) .and. (x(i)<=0.7+c)) then
                    ua(i) = 20.0*(x(i)-c)/3.0-5.0/3.0
                else
                    ua(i) = 1.0
                end if
            end do
        end if
    end subroutine

    subroutine boundvalue(u, nx)
        integer :: nx
        real*16, dimension (0:nx+1) :: u

        u(0) = u(nx-1) 
        u(nx+1) = u(2)

        un(0) = un(nx-1) 
        un(nx+1) = un(2)

    end subroutine

    subroutine output(x, un, nx, IO)
        integer :: IO
        integer :: nx
        integer :: i 
        real*16, dimension(0:nx+1) :: un
        real*16, dimension(nx) :: x 

        write(IO, *) 'variables = "x", "u"'
        write(IO, *) 'ZONE I=', nx 

        do i = 1, nx 
            write(IO, *) x(i), un(i)
        end do
    end subroutine 


end module

program convection
use convection_module
implicit none
    integer, parameter :: IO = 12
    integer :: m, nx, nt, signal_type, scheme, i, j
    real*16 :: l, c, c0, c1, cfl, pi, k, h, dt, time, phi, t, g
    real*16, allocatable, dimension(:) :: u, un, x, ua

    open(IO, file = 'input.txt')
    read(IO, *) L
    read(IO, *) m
    read(IO, *) c 
    read(IO, *) c0, c1 
    read(IO, *) nx
    read(IO, *) nt
    read(IO, *) cfl
    read(IO, *) signal_type
    read(IO, *) scheme
    close(IO)

    allocate(u(0:nx+1), un(0:nx+1), ua(0:nx+1), x(0:nx+1))

    pi=3.14159265359d0
    k = m*pi/L
    h = L/(nx-1) 
    dt = cfl*h/c 
    time = nt * dt 
    phi = -k*c/dt

    write(*,*) 'signal=', signal_type, 'scheme=', scheme
    WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
    WRITE(*,*) 'CFL=', CFL, 'dt=', dt, 'Time=', Time, 'NT=', NT
    WRITE(*,*) 'time = ', time

    x(0) = -h
    do i = 1, nx +1
        x(i) = x(i-1) + h
    end do
    
    u(:) = 0.0
    un(:) = 0.0

    t = 0.0d0
    call initvalue(un, u, ua, c, c0, c1, x, k, nx, signal_type)
    call boundvalue(u, nx)

    if (signal_type == 1) then 
        open(IO, file = 'init_1.plt')
        call output(x, un, nx, IO)
        close(IO)

    elseif (signal_type == 2) then 
        open(IO, file = 'init_2.plt')
        call output(x, un, nx, IO)
        close(IO)
    end if

    !противопоточная схема первого порядка

    if (scheme == 1) then
    g = 1.0
    write(*,*) 'g= ', g
        do j = 1, nt
            do i = 1, nx
                if (c>0) then
                    un(i) = u(i) - cfl * (u(i)-u(i-1))
                else 
                    un(i) = u(i) - cfl * (u(i+1)-u(i))
                endif
            end do
            u = un
            call boundvalue(u, un, nx)
        end do

        if (signal_type == 1) then
            open(IO, file= 'upwind_1.plt')
            call output(x, un, nx, IO)
            close(IO)
            un = c0 + c1*sin(k*(x-c*time))
            open(IO, file = 'true_1.plt')
            call output(x, un, nx, IO)
            close(IO)
            un = c0 + c1*g**nt * sin(k*(x-c*time))
            open(IO, file = 'upwind_exact_1.plt')
            call output(x, un, nx, IO)
            close(IO)

        elseif (signal_type == 2) then
            open(IO, file= 'upwind_2.plt')
            call output(x, un, nx, IO)
            close(IO)
            call initvalue(un, u, ua, c, c0, c1, x, k, nx, signal_type)
            un = ua
            open(IO, file= 'true_2.plt')
            call output(x, un, nx, IO)
            close(IO)
        endif
    endif

    !противопоточная схема второго порядка

    if (scheme == 2) then 
    g = sqrt((1-cfl+cfl*cos(k*h)*(2-cos(k*h)))**2+(cfl*sin(k*h)*(cos(k*h)-2))**2)
    write(*,*) 'g= ', g
        if (c>0) then 
            un(1) = u(1) - (cfl/2)*(3*u(1)-4*u(nx-1)+u(nx-2))
                do j = 1, nt
                    do i = 2, nx
                        un(i) = u(i)-(cfl/2)*(3*u(i)-4*u(i-1)+u(i-2))
                    end do
                    u = un
                    call boundvalue(u, un, nx)
                end do

        elseif (c<0) then
            do j = 1, nt
                do i = 1, nx-2
                    un(i) = u(i) + (cfl/2)*(3*u(i)-4*u(i+1)+u(i+2))
                end do
                un(nx-1) = u(nx) + (cfl/2)*(3*u(nx-1)-4*u(nx)+u(2))
                un(nx) = u(nx) + (cfl/2)*(3*u(nx)-4*u(2)+u(3))
                un(0) = un(nx - 1)
                un(nx+1) = un(2)
                u = un
                call boundvalue(u, un, nx)
            end do
        end if

            if (signal_type == 1) then
                open(IO, file='2nd_1.plt')
                call output(x, un, nx, IO)
                close(IO)
                un = c0 + c1*sin(k*(x-c*time))
                open(IO, file = 'true_1.plt')
                call output(x, un, nx, IO)
                close(IO)
                un = c0 + c1*g**nt * sin(k*(x-c*time))
                open(IO, file = '2nd_exact_1.plt')
                call output(x, un, nx, IO)
                close(IO)
            elseif (signal_type == 2) then
                open(IO, file= '2nd_2.plt')
                call output(x, un, nx, IO)
                close(IO)
                call initvalue(un, u, ua, c, c0, c1, x, k, nx, signal_type)
                un = ua
                open(IO, file= 'true_2.plt')
                call output(x, un, nx, IO)
                close(IO)
            endif


    end if

end program
