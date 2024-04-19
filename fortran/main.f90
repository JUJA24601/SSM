program main
    implicit none
    integer :: n = 3
    real x, y, z
    integer i
    open (17, file="test.csv", status="old")
    read (17, "()")
    do i = 1, n
        read (17, *) x, y, z
        print *, x, y, z
    end do
    close (17)
end program main