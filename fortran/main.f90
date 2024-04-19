program main
    implicit none
    integer n
    real x, y, z
    integer i

    ! csv からデータを読み込む
    open (17, file="test.csv", status="old")
    read (17, "()")
    do 
        read (17, *, end=100) x, y, z
        n = n + 1
    end do
    100 continue
    rewind(17)
    print *, "NumRec =", n
    read (17, '()')
    do i = 1, n
        read (17, *) x, y, z
        print *, x, y, z
    end do
    close (17)
end program main