module explicit_module
contains
    subroutine initvalue(y, u, un, ny, c0, c1)
    implicit none
        integer :: ny, i 
        real, dimension(1:ny) :: u, un, y
        real :: c0, c1

        do i = 1, ny
            u(i) = 0.0
            un(i) = 0.0
        end do
    end subroutine

    subroutine boundvalue(y, u, un, ny, c0, c1)
    implicit none
        integer :: ny
        real, dimension(1:ny) :: u, un, y
        real :: c0, c1 
        u(1) = 0.0
        u(ny) = 0.0
        un(1) = 0.0
        un(ny) = 0.0
    end subroutine 
end module 

program explicit_scheme
use explicit_module
implicit none
    real :: h, time, vnm, nu, w, amp, dy, dt, t_0, pi, del, c0, c1, t, tau
    real, dimension(:), allocatable :: u, un, y, u_a
    integer, parameter :: IO = 12
    integer :: ny, nt, i, m, n


    c0 = 0.0
    c1 = 0.0

    pi = 3.14159265359d0

    ny = 201
    h = 1.0
    vnm = 0.06
    amp = 1.0

    !пуазейль
    w = 0.01
    nu = 0.5
    time = 2.0
    !стокс
    ! w = 520.0
    ! nu = 20.0
    ! time = 0.05


    dy = h/(ny)
    dt = vnm*dy**2/nu
    nt = time/dt

    del = sqrt(2*nu/w)
    print*, 'delta=', del
    if (del .ge. h/2 ) then
        print*, 'low'
    else 
        print*, 'high'
    end if

    allocate(u(ny), un(ny), y(ny), u_a(ny))

    write(*,*) 'h=', h, 'dy=', dy, 'ny=', ny
    write(*,*) 'nu=', nu, 'dt=', dt, 'time =', time, 'nt=', nt

    y(1) = 0.0
    do i = 2, ny-1
        y(i) = y(i-1) + dy
    end do
    y(ny) = h

    call initvalue(y, u, un, ny, c0, c1)
    call boundvalue(y, u, un, ny, c0, c1)

    open(IO, file = 'explicit_ust.plt')
    write(IO, *) 'variables = "t", "u"'
    write(IO, *) 'ZONE I =', nt
    t = 0.0d0
    do i = 1, nt

        do m = 2, ny-1
            un(m) = u(m) + vnm*(u(m-1)-2*u(m)+u(m+1))-amp*dt*cos(w*t)
        end do
        call boundvalue(y, u, un, ny, c0, c1)
        u = un
        t = i*dt
        write(IO, *) t, u(ny/2) 
    end do
    close(IO)

    print*, ny, un(ny/2)

    open(IO, file = 'explicit.plt')
    write(IO, *) 'variables = "y", "u"'
    write(IO, *) 'ZONE I = ', ny/2
    do i = 1, ny/2
        write(IO, *) y(i), un(i)
    end do
    close(IO)

end program