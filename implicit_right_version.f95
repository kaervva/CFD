module implicit_module
contains

    !поиск коэффициентов прогонки
    subroutine find_coeff(a, b, c, ny, dt, dy, vnm)
        implicit none
        real*8, dimension(2:ny-1) :: a, b, c 
        real*8 :: dt, dy, vnm, nu
        integer ny, i 

        a(2) = 0
        b(2) = 1/dt+2*vnm/dt
        c(2) = -vnm/dt
        a(ny - 1) = -vnm/dt
        b(ny - 1) = 1/dt+2*vnm/dt
        c(ny - 1) = 0
        do i = 3, ny-2
            a(i) = -vnm/dt
            b(i) = 1/dt + 2*vnm/dt
            c(i) = a(i) 
        end do
    end subroutine find_coeff

    !процедура прогонки
    subroutine slau(ny, a, b, c, d, u)
        implicit none
        real*8, dimension(2:ny-1) :: a, b, c, d
        real*8, intent(inout) :: u(:,:)
        real*8, allocatable :: alpha(:), beta(:)
        real*8 :: k0
        integer :: i, ny
        
        allocate(alpha(3:ny-1), beta(3:ny-1))
        !прямой ход
        alpha(3) = -c(2)/b(2)
        beta(3) = d(2)/b(2)
        do i = 4, ny-1
            k0 = (b(i-1)+a(i-1)*alpha(i-1))
            alpha(i) = -c(i-1)/k0
            beta(i) = (d(i-1)-a(i-1)*beta(i-1))/k0
        end do


        u(ny-1, 2) = (d(ny-1)-a(ny-1)*beta(ny-1))/ &
                (b(ny-1)+a(ny-1)*alpha(ny-1))


        !обратный ход
        do i = ny-2, 2, -1
            u(i,2) = alpha(i+1)*u(i+1,2) + beta(i+1)
        end do
        deallocate(alpha, beta)
    end subroutine slau

    !поиск решения
    function find_result(a, b, c, d, dt, ny, nt, w, amp, vnm, IA) result(result)      
    implicit none
        real*8 :: dy, dt, t_0, h, amp, w, vnm
        real*8, dimension(2:ny-1) :: a, b, c, d
        real*8, allocatable :: u(:,:)
        real*8, allocatable :: result(:), t(:)
        integer :: i, k, ny, n, nt, m, IA
        
        allocate(u(ny,2))
        allocate(result(ny))
        allocate(t(nt))

        t_0 = 0.0
        do i = 1, nt
            t(i) = t_0
            t_0 = t_0+dt
        end do
        
        u(:,:) = 0.0

        do n = 1, nt
        
            d(2:ny-1) = u(2:ny-1,1)/dt-amp*cos(w*t(n))
            call slau(ny, a, b, c, d, u)
            write(IA,*) t(n), u(ny/2,1)
            u(:,1) = u(:,2)
        end do
        close(IA)
        result = u(:,1)

    
    end function find_result
end module implicit_module

program implicit_scheme
use implicit_module
implicit none
    real*8 :: h, w, time, vnm, dt, dy, nu, amp, dtt, t, del, tau
    real*8, allocatable, dimension (:,:) :: u
    real*8, allocatable, dimension (:) :: a, b, c, d
    real*8, allocatable, dimension (:) :: answer, y
    integer, parameter :: IO = 12, IA = 13
    integer :: ny, i, nt, m
    ny = 201
    h = 1.0
    vnm = 0.5
    amp = 1.0

    !пуазейль
    ! w = 0.01
    ! nu = 0.5
    ! time = 2.0
    !стокс
    w = 520.0
    nu = 20.0
    time = 0.03

    allocate(y(1:ny))
    allocate(answer(ny), a(2:ny-1), b(2:ny-1), c(2:ny-1), d(2:ny-1))
    allocate(u(ny, 2))

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
    print*, 'tau =', h**2/nu

    write(*,*) 'h=', h, 'dy=', dy, 'ny=', ny
    write(*,*) 'nu=', nu, 'dt=', dt, 'time =', time, 'nt=', nt

    y(1) = 0.0
    do i = 2, ny-1
        y(i) = y(i-1) + dy
    end do
    y(ny) = h

    open(IA, file='implicit_ust.plt')
    write(IA, *) 'variables = "t", "u"'
    write(IA, *) 'ZONE I =', nt

    call find_coeff(a, b, c, ny, dt, dy, vnm)

    answer = find_result(a, b, c, d, dt, ny, nt, w, amp, vnm, IA)
    print*, ny, answer(ny/2)

    open(IO, file = 'implicit.plt')
    write(IO, *) 'variables = "y", "u"'
    write(IO, *) 'ZONE I = ', ny/2
    do i = 1, ny/2
        write(IO, *) y(i), answer(i)
    end do
    close(IO)



end program