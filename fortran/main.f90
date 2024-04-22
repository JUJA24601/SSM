program main
    implicit none
    integer number_of_triangles
    double precision apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z,&
                    & a, b, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    integer i

    ! csv からデータを読み込む
    open (17, file="data.csv", status="old")
    number_of_triangles = 0
    read (17, "()")
    do 
        read (17, *, end=100) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
        a, b, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
        number_of_triangles = number_of_triangles + 1
    end do
100 continue
    rewind(17)
    print *, "number_of_triangles =", number_of_triangles
    read (17, '()')
    do i = 1, number_of_triangles
        read (17, *) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
        a, b, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    end do
    close (17)
end program main