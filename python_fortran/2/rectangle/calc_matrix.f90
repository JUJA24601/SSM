subroutine calc_matrix(rectangles_data, num_of_rectangles, P)
    implicit none
    integer i, j, num_of_rectangles
    double precision rectangles_data(num_of_rectangles,4,3), P(num_of_rectangles,num_of_rectangles), eps0, pi, coord(4,3),&
                    & a, b, p0(1,3), n1(1,3), n2(1,3), n3(1,3), E0(1,3)
    intent(in) rectangles_data, num_of_rectangles
    intent(out) P

    ! 4角形の型
    type Rectangle
        double precision coord(4,3), a, b, p0(1,3), n1(1,3), n2(1,3), n3(1,3)
    end type Rectangle
    type(Rectangle) rectangles(num_of_rectangles)  ! 4角形の構造体を収める配列

    eps0 = 8.8541878128d-12     ! 誘電率
    pi = acos(-1.0d0)           ! 円周率
    print *, "----- fortran (calc_matrix) -----"
    print *, "num_of_rectangles : ", num_of_rectangles
    ! print *, "E0 : ", E0
    
    ! 4角形の型を収める配列の初期化
    do i = 1, num_of_rectangles
        coord = rectangles_data(i,:,:)
        a = distance(coord(1,:), coord(2,:))
        b = distance(coord(1,:), coord(4,:))
        p0(1,:) = coord(1,:)
        n1(1,:) = (coord(2,:) - coord(1,:))/a
        n2(1,:) = (coord(4,:) - coord(1,:))/b
        n3 = cross_product(n1, n2)
        rectangles(i)%coord = coord
        rectangles(i)%a = a
        rectangles(i)%b = b
        rectangles(i)%p0 = p0
        rectangles(i)%n1 = n1
        rectangles(i)%n2 = n2
        rectangles(i)%n3 = n3
    end do


    ! ここから行列を計算する
    do i = 1, num_of_rectangles
        do j = 1, num_of_rectangles
            P(j,i) = p_calc(rectangles(i), point(rectangles(j)))
        end do
    end do

    print *, "calc_matrix -- done"

contains

    ! function potential_elem(tri, r, sigma)
    !     !
    !     !
    !     ! 三角形要素triangleによってできる任意の点rでの電位を求める．
    !     ! 
    !     ! Parameters
    !     ! ---------
    !     ! tri      : [3,3]の実数配列でそれぞれの行に頂点が保存される．
    !     !            三角形の座標
    !     ! r        : [1,3]の実数配列
    !     !            任意の点
    !     ! sigma    : 実数
    !     !            三角形の電荷密度
    !     !
    !     ! Returns
    !     ! ---------
    !     ! potential_elem : 実数
    !     !                  三角形要素によって点rで生じる電位
    !     ! 
    !     !
    !     implicit none
    !     double precision potential_elem, r(1,3), sigma
    !     type(Rectangle) tri
    !     potential_elem = (p_calc(tri, r)*sigma)/(4.0d0*pi*eps0)
    ! end function potential_elem

    function point(rect)
        implicit none
        double precision point(1,3)
        type(Rectangle) rect
        a  = rect%a
        b  = rect%b
        p0 = rect%p0
        n1 = rect%n1
        n2 = rect%n2
        n3 = rect%n3
        point = p0 + 0.5*(a*n1 + b*n2)
    end function point


    function distance(a, b)     ! OK
        !
        !
        ! 2点間 a,b の距離を求める．
        !
        ! Parameters
        ! ----------
        ! a : [1,3]の実数配列
        !     aの座標
        ! b : [1,3]の実数配列
        !     bの座標
        !
        ! Returns
        ! ----------
        ! distance : 実数
        !            2点間の距離
        !
        !
        implicit none
        double precision a(3), b(3), distance
        distance = sqrt(sum((a-b)*(a-b)))
    end function distance


    function p_calc(rect, r)
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
        double precision p_calc, r(1,3), a, b, p0(1,3), n1(1,3), n2(1,3), n3(1,3), xmin, xmax,&
                        & ymin, ymax, p(1,3), uP, vP, w
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

        p_calc = Integral_ln(xmax, ymax, w) - Integral_ln(xmin, ymax, w) - Integral_ln(xmax, ymin, w) + Integral_ln(xmin, ymin, w)
        
    end function p_calc

    function Integral_ln(x, y, w)
        implicit none
        double precision x, y, w, Integral_ln, r, r0, xa, c1, c2, c3
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
        Integral_ln = c1 + c2 - xa + c3
        if (x<0.0d0) then
            Integral_ln = -Integral_ln
        end if
    end function Integral_ln
    

    function cross_product(a,b)
        implicit none
        double precision a(3), b(3), cross_product(1,3)
        cross_product(1,1) = a(2)*b(3) - a(3)*b(2)
        cross_product(1,2) = a(3)*b(1) - a(1)*b(3)
        cross_product(1,3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    function area(rect)
        implicit none
        double precision area
        type(Rectangle) rect
        a  = rect%a
        b  = rect%b
        area = a*b
    end function area

end subroutine calc_matrix