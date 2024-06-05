module burgers_module
contains

    subroutine initvalue(un, u, x, k, nx, signal_type)
        integer :: nx, signal_type, i 
        real*16 :: k, c0, c1
        real*16, dimension(nx) :: x
        real*16, dimension(0:nx+1) :: un, u

        if (signal_type == 1) then
            un = sin(k*x)
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
        end if
    end subroutine

    subroutine boundvalue(u, un, nx)
        integer :: nx
        real*16, dimension (0:nx+1) :: u, un

        u(0) = u(nx-1) 
        u(nx+1) = u(2)

        un(0) = un(nx-1) 
        un(nx+1) = un(2)

    end subroutine

    subroutine output(x, u, nx, IO)
        integer:: IO, nx, i 
        real*16, dimension(0:nx+1) :: u, x

        write(IO, *) 'variables = "x", "u"'
        write(IO, *) 'ZONE I=', nx 

        do i = 1, nx 
            write(IO, *) x(i), u(i)
        end do
    end subroutine

end module

program burgers 
use burgers_module
implicit none
    integer, parameter :: IO = 12
    integer ::  m, nx, nt, signal_type, scheme, i, j
    real*16, allocatable :: u(:), un(:), x(:)
    real*16 :: l, cfl, pi, k, h, dt, time, t, c

    open(IO, file = 'input.txt')
    read(IO, *) L
    read(IO, *) m
    read(IO, *) c
    read(IO, *) nx
    read(IO, *) nt
    read(IO, *) cfl
    read(IO, *) signal_type
    read(IO, *) scheme
    close(IO)

    allocate(u(0:nx+1), un(0:nx+1), x(0:nx+1))

    pi=3.14159265359d0
    k = m*pi/L
    h = L/(nx-1) 
    dt = cfl*h
    time = nt * dt 

    WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
    WRITE(*,*) 'CFL=', CFL, 'dt=', dt, 'Time=', Time, 'NT=', NT

    x(0) = -h
    do i = 1, nx+1
        x(i) = x(i-1)+h
    end do

    u = 0.0
    un = 0.0

    t = 0.0d0 
    call initvalue(un, u, x, k, nx, signal_type)
    call boundvalue(u, un, nx)
    if (signal_type == 1) then 
        open(IO, file = 'init_1.plt')
        call output(x, u, nx, IO)
        close(IO)
    else 
        open(IO, file = 'init_2.plt')
        call output(x, u, nx, IO)
        close(IO) 
    end if 

    if (scheme == 1) then 
        do j = 1, nt 
            do i = 2, nx-1 
                if (u(i)>=0) then 
                    un(i) = u(i) - (cfl/2)*(u(i)**2-u(i-1)**2) 
                else 
                    un(i) = u(i) - (cfl/2)*(u(i+1)**2-u(i)**2) 
                end if 
            end do 
            u = un 
            call boundvalue(u, un, nx)
        end do 

        if (signal_type == 1) then 
            open(IO, file= '1st_1.plt')
            call output(x, u, nx, IO)
            close(IO)
        elseif (signal_type == 2) then
            open(IO, file= '1st_2.plt')
            call output(x, u, nx, IO)
            close(IO)
        endif
    elseif (scheme == 2) then 
        if (c>0) then
                do j = 1, nt
                un(1) = u(1) - (cfl/4)*(3*u(1)**2-4*u(nx-1)**2+u(nx-2)**2)
                    do i = 2, nx
                        un(i) = u(i)-(cfl/4)*(3*u(i)**2-4*u(i-1)**2+u(i-2)**2)
                    end do
                    u = un
                    call boundvalue(u, un, nx)
                end do
        elseif (c<0) then
            do j = 1, nt
                do i = 1, nx-2
                    un(i) = u(i) + (cfl/4)*(3*u(i)**2-4*u(i+1)**2+u(i+2)**2)
                end do
                un(nx-1) = u(nx) + (cfl/4)*(3*u(nx-1)**2-4*u(nx)**2+u(2)**2)
                un(nx) = u(nx) + (cfl/4)*(3*u(nx)**2-4*u(2)**2+u(3)**2)
                u = un
                call boundvalue(u, un, nx)
            end do
        end if
        if (signal_type == 1) then
            open(IO, file='2nd_1.plt')
            call output(x, un, nx, IO)
            close(IO)
        elseif (signal_type == 2) then
            open(IO, file= '2nd_2.plt')
            call output(x, un, nx, IO)
            close(IO)
        endif
    endif
end program
