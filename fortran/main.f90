program main
    implicit none
    integer i, j, num_of_rectangles, num_of_triangles, info
    double precision apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z,&
                    apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    double precision potential_conductor, Q, C, pi, eps0

    ! type of rectangle
    type Rectangle
        double precision coord(4,3), a, b, p0(3), n1(3), n2(3), n3(3)
    end type Rectangle

    ! type of triangle
    type Triangle
        double precision coord(3,3), a, b, p0(3), n1(3), n2(3), n3(3)
    end type Triangle

    type(Rectangle), allocatable :: rectangles(:)
    double precision, allocatable :: P(:,:)
    double precision, allocatable :: ipiv(:)
    double precision, allocatable :: sigmas(:)

    potential_conductor = 10    ! potential of conductor
    pi = acos(-1.0d0)
    eps0 = 8.8541878128d-12     ! permittivity

    ! ------------ read .csv file : Rectangles ------------
    ! open (17, file="data.csv", status="old")
    ! num_of_rectangles = 0
    ! read (17, "()")
    ! do
    !     read (17, *, end=100) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
    !     apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    !     num_of_rectangles = num_of_rectangles + 1
    ! end do
    ! 100 continue
    ! rewind(17)
    ! allocate(rectangles(num_of_rectangles))
    ! allocate(P(num_of_rectangles, num_of_rectangles))
    ! allocate(ipiv(num_of_rectangles))
    ! allocate(sigmas(num_of_rectangles))
    ! print *, "num_of_rectangles =", num_of_rectangles
    ! read (17, '()')
    ! do i = 1, num_of_rectangles
    !     read (17, *) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
    !     apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    !     rectangles(i)%coord(1,:) = [apex1_x, apex1_y, apex1_z]
    !     rectangles(i)%coord(2,:) = [apex2_x, apex2_y, apex2_z]
    !     rectangles(i)%coord(3,:) = [apex3_x, apex3_y, apex3_z]
    !     rectangles(i)%coord(4,:) = [apex4_x, apex4_y, apex4_z]
    !     rectangles(i)%a          = a_temp
    !     rectangles(i)%b          = b_temp
    !     rectangles(i)%p0(:) = [p0_x, p0_y, p0_z]
    !     rectangles(i)%n1(:) = [n1_x, n1_y, n1_z]
    !     rectangles(i)%n2(:) = [n2_x, n2_y, n2_z]
    !     rectangles(i)%n3(:) = [n3_x, n3_y, n3_z]
    ! end do
    ! close (17)
    
    ! ------------ Triangles ------------
    open (17, file="data.csv", status="old")
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
    allocate(P(num_of_triangles, num_of_triangles))
    allocate(ipiv(num_of_triangles))
    allocate(sigmas(num_of_triangles))
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
        triangles(i)%p0(:) = [p0_x, p0_y, p0_z]
        triangles(i)%n1(:) = [n1_x, n1_y, n1_z]
        triangles(i)%n2(:) = [n2_x, n2_y, n2_z]
        triangles(i)%n3(:) = [n3_x, n3_y, n3_z]
    end do
    close (17)

    ! calculate coefficient matrix ------------ Rectangle ------------
    ! do i = 1, num_of_rectangles
    !     do j = 1, num_of_rectangles
    !         P(j,i) = p_calc_rect(rectangles(i), centroid_rectangle(rectangles(j)))
    !     end do
    ! end do

    ! ------------ Triangle ------------
    do i = 1, num_of_triangles
        do j = 1, num_of_triangles
            P(j,i) = p_calc_tri(triangles(i), centroid_triangle(triangles(j)))
        end do
    end do

    ! RHS ------------ Rectangle ------------
    ! do i = 1, num_of_rectangles
    !     sigmas(i) = potential_conductor
    ! end do

    ! ------------ Triangle ------------
    do i = 1, num_of_triangles
        sigmas(i) = potential_conductor
    end do

    ! solve linear system ------------ Rectangle ------------
    ! call dgesv(num_of_rectangles, 1, P, num_of_rectangles, ipiv, sigmas, num_of_rectangles, info)
    ! sigmas = 4*pi*eps0*sigmas
    ! print *, "sigmas : ", sigmas(1)

    ! ------------ Triangle ------------
    call dgesv(num_of_triangles, 1, P, num_of_triangles, ipiv, sigmas, num_of_triangles, info)
    sigmas = 4*pi*eps0*sigmas
    print *, "sigmas : ", sigmas(1)

    ! calculate total charge ------------ Rectangle ------------
    ! Q = 0
    ! do i = 1, num_of_rectangles
    !     Q = Q + sigmas(i)*area_rect(rectangles(i))
    ! end do
    ! print *, "total charge : ", Q

    ! ------------ Triangle ------------
    Q = 0
    do i = 1, num_of_triangles
        Q = Q + sigmas(i)*area_tri(triangles(i))
    end do
    print *, "total charge : ", Q

    ! calculate capacitance
    C = Q/potential_conductor
    print *, "Capacitance : ", C
    print *, "Capacitance/(4*pi*eps0)", C/(4*pi*eps0)


contains

    function centroid_rectangle(rect) result(centroid)
        !
        !
        ! 四角形の中心を求める
        ! rect      : 四角形の型
        ! centorid  : 四角形の中心
        !
        !
        implicit none
        type(Rectangle), intent(in) :: rect
        double precision centroid(3), a, b, p0(3), n1(3), n2(3)
        a  = rect%a
        b  = rect%b
        p0 = rect%p0
        n1 = rect%n1
        n2 = rect%n2
        centroid = p0 + 0.5*(a*n1 + b*n2)
    end function

    function centroid_triangle(triangle) result(centroid)
        !
        !
        ! 三角形の内心を求める．
        !
        ! Parameters
        ! -----------
        ! triangle : [3,3]の実数配列でそれぞれの行に頂点が保存される
        !            三角形の座標
        !
        ! Returns
        ! -----------
        ! incenter : [1,3]の実数配列
        !            三角形の内心
        !
        !
        implicit none
        type(Triangle), intent(in) :: triangle
        double precision triangle(3,3), centroid(3), a(3), b(3), c(3), ab, bc, ca
        a = triangle(1,:); b = triangle(2,:); c = triangle(3,:)
        ab = distance(a,b); bc = distance(b,c); ca = distance(c,a)
        centroid = (bc*a+ca*b+ab*c)/(ab+bc+ca)
    end function


    function distance(a, b) result(dist)
        !
        !
        ! 2点間 a,b の距離を求める．
        !
        !
        implicit none
        double precision a(3), b(3), dist
        dist = sqrt(sum((a-b)*(a-b)))
    end function


    function p_calc_rect(rect, r) result(ans)
        ! 
        !
        ! 電位係数を求める．(4角形がrに作る)
        !
        ! Parameters
        ! -----------
        ! rectangle : 4角形の構造体
        ! r   : 任意の点
        !
        ! Returns
        ! -----------
        ! p_calc : 実数
        !          電位係数
        !
        !
        implicit none
        double precision ans, r(3), a, b, p0(3), n1(3), n2(3), n3(3), xmin, xmax,&
                        & ymin, ymax, p(3), uP, vP, w
        type(Rectangle) rect

        a  = rect%a
        b  = rect%b
        p0 = rect%p0
        n1 = rect%n1
        n2 = rect%n2
        n3 = rect%n3

        p  = r - p0
        uP = sum(p*n1)
        vP = sum(p*n2)
        w  = sum(p*n3)

        xmin = -uP
        xmax = -uP + a
        ymin = -vP
        ymax = -vP + b

        ans = Integral_ln(xmax, ymax, w) - Integral_ln(xmin, ymax, w) - Integral_ln(xmax, ymin, w) + Integral_ln(xmin, ymin, w)
        
    end function

    function Integral_ln(x, y, w) result(ans)
        implicit none
        double precision x, y, w, ans, r, r0, xa, c1, c2, c3
        r  = sqrt(abs(x*x + y*y + w*w))
        r0 = sqrt(abs(y*y + w*w))
        xa = abs(x)
        if (xa < 1.0d-10) then
            c1 = 0.0d0
        else
            c1 = xa*log(abs(y + r) + 1.0d-12)
        end if
        if (abs(y) < 1.0d-12) then
            c2 = 0.0d0
        else
            c2 = y*log(abs((xa + r)/r0) + 1.0d-12)
        end if
        if (abs(w) < 1.0d-12) then
            c3 = 0.0d0
        else
            c3 = w*(atan(xa/w) + atan(y*w/(x*x + w*w + xa*r)) - atan(y/w))
        end if
        ans = c1 + c2 - xa + c3
        if (x<0.0d0) then
            ans = -ans
        end if
    end function
    
    function p_calc_tri(tri, r) result(ans) ! ここの変更
        ! 
        !
        ! 電位係数を求める．(三角形がrに作る)
        !
        ! Parameters
        ! -----------
        ! tri : 三角形の構造体
        ! r   : 任意の点
        !
        ! Returns
        ! -----------
        ! p_calc : 実数
        !          電位係数
        !
        !
        implicit none
        double precision p_calc, r(3), a, b, p0(3), n1(3), n2(3), n3(3), x1, x2, x3,&
                        & y1, y2, z, a1, a2, b1, b2, u1, u2, p(1,3), n1dotn2, N2prime(3), n2dotn2prime
        type(Triangle) tri

        ! 三角形の情報を抜き出す
        a  = tri%a
        b  = tri%b
        p0 = tri%p0
        n1 = tri%n1
        n2 = tri%n2
        n3 = tri%n3
        
        p = p0 - r

        n1dotn2 = sum(n1*n2)
        N2prime = (n2 - (n1dotn2 * n1))/sqrt(sum((n2 - (n1dotn2 * n1))*(n2 - (n1dotn2 * n1))))
        n2dotn2prime = sum(n2*N2prime)
        x1 = sum(p*n1)
        y1 = sum(p*N2prime)
        z = -sum(p*n3)

        x2 = x1 + a

        x3 = x1 + b * n1dotn2

        if (z < 0.0d0) then
            y1 = -y1
            y2 = y1 - b * n2dotn2prime
        else
            y2 = y1 + b * n2dotn2prime
        end if

        z = abs(z)

        if (z > 1.0d-14) then
            u1 = y1 / z
            if (u1 == 0.0d0) then
                u1 = 1.0d-18
            end if
            u2 = y2 / z
            if (u2 == 0.0d0) then
                u2 = 1.0d-18
            end if
            a1 = (x1 * y2 - x3 * y1) / (z * (y2 - y1))
            a2 = (x2 * y2 - x3 * y1) / (z * (y2 - y1))
        else
            u1 = y1
            if (u1 == 0.0d0) then
                u1 = 1.0d-18
            end if
            u2 = y2
            if (u2 == 0.0d0) then
                u2 = 1.0d-18
            end if
            a1 = (x1 * y2 - x3 * y1) / (y2 - y1)
            a2 = (x2 * y2 - x3 * y1) / (y2 - y1)
        end if

        b1 = (x3 - x1) / (y2 - y1)
        b2 = (x3 - x2) / (y2 - y1)

        if (z > 1.0d-14) then
            if (abs(b1) < 1.0d-13) then
                p_calc = z * (I1(a2, b2, u1, u2) - I2(a1, u1, u2))
            else if (abs(b2) < 1.0d-13) then
                p_calc = z * (I2(a2, u1, u2) - I1(a1, b1, u1, u2))
            else
                p_calc = z * (I1(a2, b2, u1, u2) - I1(a1, b1, u1, u2))
            end if
        else
            p_calc = (p_calc_noZ(a2, b2, a1, b1, u2) - p_calc_noZ(a2, b2, a1, b1, u1))
        end if

        p_calc = abs(p_calc)
        
    end function p_calc

    ! function cross_product(a,b)
    !     implicit none
    !     double precision a(3), b(3), cross_product(1,3)
    !     cross_product(1,1) = a(2)*b(3) - a(3)*b(2)
    !     cross_product(1,2) = a(3)*b(1) - a(1)*b(3)
    !     cross_product(1,3) = a(1)*b(2) - a(2)*b(1)
    ! end function cross_product

    function area_rect(rect) result(ans)
        implicit none
        double precision a, b, ans
        type(Rectangle) rect
        a  = rect%a
        b  = rect%b
        ans = a*b
    end function

    function area_tri(triangle) result(area)
        implicit none
        double precision triangle(3,3), area, a(3), b(3), o(3)
        o = [0.0d0, 0.0d0, 0.0d0]
        a = triangle(2,:) - triangle(1,:)
        b = triangle(3,:) - triangle(1,:)
        area = 0.5*sqrt((distance(a,o)**2)*(distance(b,o)**2) - sum(a*b)**2)
    end function

end program main