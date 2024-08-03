program main
    implicit none
    integer i, j, num_of_rectangles, info
    double precision apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z,&
                    apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
    double precision potential_conductor, Q, C, pi, eps0

    ! derive type of rectangle
    type Rectangle
        double precision coord(4,3), a, b, p0(3), n1(3), n2(3), n3(3)
    end type Rectangle

    type(Rectangle), allocatable :: rectangles(:)
    double precision, allocatable :: P(:,:)
    double precision, allocatable :: ipiv(:)
    double precision, allocatable :: sigmas(:)

    potential_conductor = 10    ! potential of conductor
    pi = acos(-1.0d0)
    eps0 = 8.8541878128d-12     ! permittivity

    ! read .csv file
    open (17, file="data.csv", status="old")
    num_of_rectangles = 0
    read (17, "()")
    do
        read (17, *, end=100) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
        apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
        num_of_rectangles = num_of_rectangles + 1
    end do
    100 continue
    rewind(17)
    allocate(rectangles(num_of_rectangles))
    allocate(P(num_of_rectangles, num_of_rectangles))
    allocate(ipiv(num_of_rectangles))
    allocate(sigmas(num_of_rectangles))
    print *, "num_of_rectangles =", num_of_rectangles
    read (17, '()')
    do i = 1, num_of_rectangles
        read (17, *) apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, &
        apex4_x, apex4_y, apex4_z, a_temp, b_temp, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z
        rectangles(i)%coord(1,:) = [apex1_x, apex1_y, apex1_z]
        rectangles(i)%coord(2,:) = [apex2_x, apex2_y, apex2_z]
        rectangles(i)%coord(3,:) = [apex3_x, apex3_y, apex3_z]
        rectangles(i)%coord(4,:) = [apex4_x, apex4_y, apex4_z]
        rectangles(i)%a          = a_temp
        rectangles(i)%b          = b_temp
        rectangles(i)%p0(:) = [p0_x, p0_y, p0_z]
        rectangles(i)%n1(:) = [n1_x, n1_y, n1_z]
        rectangles(i)%n2(:) = [n2_x, n2_y, n2_z]
        rectangles(i)%n3(:) = [n3_x, n3_y, n3_z]
    end do
    close (17)
    

    ! ! 確認用
    ! print *, "centroid_rectangle : "
    ! print *, centroid_rectangle(rectangles(1))
    ! print *, "distance : "
    ! print *, distance(rectangles(1)%coord(1,:), rectangles(1)%coord(3,:))
    ! print *, "p_calc_rect : "
    ! print *, p_calc_rect(rectangles(1), centroid_rectangle(rectangles(1)))
    ! print *, "area_rect : "
    ! print *, area_rect(rectangles(1))

    ! calculate coefficient matrix
    do i = 1, num_of_rectangles
        do j = 1, num_of_rectangles
            P(j,i) = p_calc_rect(rectangles(i), centroid_rectangle(rectangles(j)))
        end do
    end do

    ! RHS
    do i = 1, num_of_rectangles
        sigmas(i) = potential_conductor
    end do

    ! solve linear system
    call dgesv(num_of_rectangles, 1, P, num_of_rectangles, ipiv, sigmas, num_of_rectangles, info)
    sigmas = 4*pi*eps0*sigmas
    print *, "sigmas : ", sigmas(1)

    ! calculate total charge
    Q = 0
    do i = 1, num_of_rectangles
        Q = Q + sigmas(i)*area_rect(rectangles(i))
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

end program main