module TRIANGLE_mod
    implicit none
    type Triangle
        double precision coord(3,3), a, b, p0(1,3), n1(1,3), n2(1,3), n3(1,3)
    end type Triangle
end module TRIANGLE_mod

program main
    use TRIANGLE_mod
    implicit none
    integer num_of_triangles
    double precision apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z,&
                    & a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    integer i
    type(Triangle), allocatable :: triangles(:)

    ! csv からデータを読み込む
    open (17, file="test.csv", status="old")
    num_of_triangles = 0
    read (17, "()")
    do 
        read (17, *, end=100) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
        a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
        num_of_triangles = num_of_triangles + 1
    end do

100 continue
    rewind(17)
    allocate(triangles(num_of_triangles))
    print *, "num_of_triangles =", num_of_triangles
    read (17, '()')
    do i = 1, num_of_triangles
        read (17, *) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
        a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
        triangles(i)%coord(1,:) = [apex1_x, apex1_y, apex1_z]
        triangles(i)%coord(2,:) = [apex2_x, apex2_y, apex2_z]
        triangles(i)%coord(3,:) = [apex3_x, apex3_y, apex3_z]
        triangles(i)%a          = a_temp
        triangles(i)%b          = b_temp
        triangles(i)%p0(1,:) = [p0_x, p0_y, p0_z]
        triangles(i)%n1(1,:) = [n1_x, n1_y, n1_z]
        triangles(i)%n2(1,:) = [n2_x, n2_y, n2_z]
        triangles(i)%n3(1,:) = [n3_x, n3_y, n3_z]
    end do

    print *, "triangles(1) =", triangles(1)%coord
    print *, "triangles(2) =", triangles(2)%coord


    close (17)
    deallocate(triangles)
end program main