program Fractal

    use omp_lib
    use plplot
    implicit none

    real(8), parameter :: CMAP_SCALE = 0.125_8
    real(8), parameter :: CMAP_OFFSET = 0

    integer :: width, height
    real(8) :: zoom, bailout
    integer :: max_iter
    real(10) :: x, y, dx, dy
    integer :: num_threads
    real(8) :: start, finish, time
    real(8), allocatable :: image(:, :)

    print *, 'Geometry (width, height):'
    read *, width, height
    print *, 'Zoom level (in powers of 2):'
    read *, zoom
    print *, 'Bailout radius:'
    read *, bailout
    print *, 'Maximum number of iterations:'
    read *, max_iter

    if (zoom > 0) then
        ! Interesting region
        x = -0.74364388703715870475_10
        y =  0.131825904205311970493_10
    else
        x = -0.5_10
        y = 0
    end if
    dx = 2 * 2.0_10**(-zoom)
    dy = dx * height / width

    allocate (image(width, height))

!    call omp_set_num_threads(1)
    num_threads = omp_get_max_threads()

    call cpu_time(start)
    call mandelbrot(x - dx, x + dx, y - dy, y + dy, bailout, max_iter, image)
    call cpu_time(finish)

    time = (finish - start) / num_threads
    print '("Elapsed time: ", f6.3, " s (", i1, " threads)")', time, num_threads

    call filter(image, CMAP_SCALE, CMAP_OFFSET)
    call plot_fractal(image)

    deallocate (image)

contains

    subroutine mandelbrot(min_x, max_x, min_y, max_y, bailout, max_iter, image)
        real(10), intent(in) :: min_x, max_x, min_y, max_y
        real(8), intent(in) :: bailout
        integer, intent(in) :: max_iter
        real(8), intent(out) :: image(:, :)

        integer :: num_x, num_y
        integer :: i, j
        real(10) :: x, y

        num_x = size(image, 1)
        num_y = size(image, 2)
        image = 0
        !$omp parallel do schedule(static, 1) private(i, j, x, y)
        do i = 1, num_x
            x = min_x + (i - 1) * (max_x - min_x) / (num_x - 1)
            do j = 1, num_y
                y = min_y + (j - 1) * (max_y - min_y) / (num_y - 1)
                image(i, j) = escape_time(x, y, bailout, max_iter)
            end do
        end do
        !$omp end parallel do

    end subroutine

    real(8) function escape_time(x, y, bailout, max_iter) result(n)
        real(10), intent(in) :: x, y
        real(8), intent(in) :: bailout
        integer, intent(in) :: max_iter

        integer  :: i
        complex(10) :: c, z

        n = 0

        if (is_in_bulb(x, y) .or. is_in_cardioid(x, y)) return

        c = cmplx(x, y)
        z = c
        do i = 0, max_iter - 1
            if (abs(z) > bailout) then
                n = dble(i - log(log(abs(z)) / log(bailout)) / log(2.0_10))
                return
            end if
            z = z**2 + c
        end do
    end function

    logical function is_in_cardioid(x, y)
        real(10), intent(in) :: x, y

        real(10) :: q

        q = (x - 0.25_10) * (x - 0.25_10) + y * y
        is_in_cardioid = q * (q + (x - 0.25_10)) < 0.25_10 * y * y
    end function

    logical function is_in_bulb(x, y)
        real(10), intent(in) :: x, y

        is_in_bulb = (x + 1) * (x + 1) + y * y < 1.0_10 / 16
    end function

    subroutine filter(image, cmap_scale, cmap_offset)
        real(8), intent(inout) :: image(:, :)
        real(8), intent(in) :: cmap_scale, cmap_offset

        image = modulo(sqrt(image) * cmap_scale + cmap_offset, 1.0_8)
    end subroutine

    subroutine plot_fractal(image)
        real(8), intent(in) :: image(:, :)

        character(13) :: geom
        real(8) :: aspect

        call plparseopts(PL_PARSE_FULL)

        write (geom, '(i6, "x", i6)') size(image, 1), size(image, 2)
        call plsetopt('geometry', geom)

        call plinit()

        aspect = 1.0_8 * size(image, 1) / size(image, 2)
        call plenv(0.0_8, aspect, 0.0_8, 1.0_8, 1, -2)
        call pllab('Re[c]', 'Im[c]', 'Mandelbrot fractal')

        call set_cmap(4096)

        call plimage(image, 0.0_8, aspect, 0.0_8, 1.0_8, 0.0_8, 1.0_8, 0.0_8, aspect, 0.0_8, 1.0_8)

        call plend()
    end subroutine

    subroutine set_cmap(ncol)
        integer, intent(in) :: ncol

        real(8), parameter :: X(7) = [ 0.0_8,  3.0_8,  6.0_8,  9.0_8, 10.0_8, 11.0_8, 12.0_8] / 12
        real(8), parameter :: R(7) = [ 0.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8,  1.0_8,  0.0_8]
        real(8), parameter :: G(7) = [ 0.0_8,  0.0_8,  1.0_8,  1.0_8,  1.0_8,  0.0_8,  0.0_8]
        real(8), parameter :: B(7) = [ 0.0_8,  0.0_8,  0.0_8,  1.0_8,  0.0_8,  0.0_8,  0.0_8]

        call plscmap1n(ncol)
        call plscmap1l(.true., X, R, G, B)
    end subroutine

end program

