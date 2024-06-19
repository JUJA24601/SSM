subroutine calc_efield(triangles_data, sigmas, E0, num_of_triangles, pi)
    implicit none
    integer i, j, k, num_of_triangles
    double precision triangles_data(num_of_triangles,3,3), eps0, pi, coord(3,3), side(3),&
                    & a, b, p0(1,3), n1(1,3), n2(1,3), n3(1,3), E0(1,3), r(1,3), sigmas(num_of_triangles),&
                    & f(1,3)
    intent(in) triangles_data, sigmas, E0, num_of_triangles
    intent(out) pi

    ! 三角形の型
    type Triangle
        double precision coord(3,3), a, b, p0(1,3), n1(1,3), n2(1,3), n3(1,3)
    end type Triangle
    type(Triangle) triangles(num_of_triangles)  ! 三角形の構造体を収める配列

    eps0 = 8.8541878128d-12     ! 誘電率
    pi = acos(-1.0d0)           ! 円周率
    print *, "----- fortran (calc_Efield) -----"
    print *, "num_of_triangles : ", num_of_triangles
    print *, "E0 : ", E0
    
    ! 三角形の型を収める配列の初期化
    do i = 1, num_of_triangles
        coord = triangles_data(i,:,:)
        side = [distance(coord(1,:),coord(2,:)), distance(coord(2,:), coord(3,:)), distance(coord(3,:), coord(1,:))]
        a = MAXVAL(side)
        p0 = triangles_data(i, MAXLOC(side),:)
        b = distance(coord(mod((MAXLOC(side)+1),3)+1,:), p0)
        n1 = (coord(mod(MAXLOC(side),3)+1,:) - p0)/distance(coord(mod(MAXLOC(side),3)+1,:),p0)
        n2 = (coord(mod((MAXLOC(side)+1),3)+1,:) - p0)/distance(coord(mod((MAXLOC(side)+1),3)+1,:),p0)
        n3 = cross_product(n1,n2)/sqrt(sum((cross_product(n1,n2))*(cross_product(n1,n2))))
        triangles(i)%coord = coord
        triangles(i)%a = a
        triangles(i)%b = b
        triangles(i)%p0 = p0
        triangles(i)%n1 = n1
        triangles(i)%n2 = n2
        triangles(i)%n3 = n3
    end do


    ! ここから電場を求める
    ! do i = 1, 1
    !     do j = 1, 1
    !         do k = 1, 1
    !             ! r(1,:) = centroid(triangles(1)%coord)
    !             ! r(1,3) = 10.0d0
    !             ! ! print *, r
    !             ! f = field_elem(triangles(1), r, 1.0d-8)       ! 外部電場がある場合にはここに外部電場を追加(今後)
    !             ! print *, "f : ", f

    !             ! r(1,:) = [sqrt(10.0d0), 0.0d0, 0.0d0]
    !             r(1,:) = [0.0d0, 0.0d0, 0.0d0]
    !             f = field(triangles, r, sigmas, num_of_triangles)       ! 外部電場がある場合にはここに外部電場を追加(今後)
    !             print *, "f : ", f
    !         end do
    !     end do
    ! end do
    
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
    !     type(Triangle) tri
    !     potential_elem = (p_calc(tri, r)*sigma)/(4.0d0*pi*eps0)
    ! end function potential_elem


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


    function area(triangle)     ! OK
        !
        ! 三角形の面積を求める．
        !
        ! Parameters
        ! -----------
        ! triangle : [3,3]の実数配列でそれぞれの行に頂点が保存される．
        !            三角形の座標
        !
        ! Returns
        ! -----------
        ! area  : 実数
        !         三角形の面積
        !
        !
        implicit none
        double precision triangle(3,3), area, a(3), b(3), o(3)
        o = [0.0d0, 0.0d0, 0.0d0]
        a = triangle(2,:) - triangle(1,:)
        b = triangle(3,:) - triangle(1,:)
        area = 0.5*sqrt((distance(a,o)**2)*(distance(b,o)**2) - sum(a*b)**2)
        
    end function area

    function incenter(triangle)     ! OK
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
        double precision triangle(3,3), incenter(3), a(3), b(3), c(3), ab, bc, ca
        a = triangle(1,:); b = triangle(2,:); c = triangle(3,:)
        ab = distance(a,b); bc = distance(b,c); ca = distance(c,a)
        incenter = (bc*a+ca*b+ab*c)/(ab+bc+ca)
    end function incenter

    function field(triangles, r, sigmas, num_of_triangles)
        implicit none
        integer num_of_triangles, i
        double precision r(3), sigmas(num_of_triangles), field(1,3)
        type(Triangle) triangles(num_of_triangles)
        field(1,:) = [0.0d0, 0.0d0, 0.0d0]
        do i = 1, num_of_triangles
            field = field + field_elem(triangles(i), r, sigmas(i))
        end do
    end function field


    function field_elem(tri, r, sigma)
        ! 
        !
        ! 電場を求める．(三角形がrに作る)
        !
        ! Parameters
        ! -----------
        ! tri : 三角形の構造体
        ! r   : 任意の点
        !
        ! Returns
        ! -----------
        ! field : 電場
        !        
        !
        !
        implicit none
        double precision field_elem(1,3), r(1,3), a, b, p0(1,3), n1(1,3), n2(1,3), n3(1,3), x1, x2, x3,&
                        & y1, y2, z, a1, a2, b1, b2, u1, u2, p(1,3), n1dotn2, N2prime(1,3), n2dotn2prime,&
                        & z_sign, local_field(1,3), dist, prefac, sigma
        type(Triangle) tri

        ! 三角形の情報を抜き出す
        a  = tri%a
        b  = tri%b
        p0 = tri%p0
        n1 = tri%n1
        n2 = tri%n2
        n3 = tri%n3
        
        dist = distance(centroid(tri%coord), r)
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
            z_sign = -1
        else
            y2 = y1 + b * n2dotn2prime
            z_sign = 1
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

        prefac = 1.0d0/(4*pi*eps0)

        if (z > 1.0d-14) then
            local_field(1,1) = z_sign * prefac * local_Ex(a1, a2, b1, b2, u1, u2)
            local_field(1,2) = prefac * local_Ey(a1, a2, b1, b2, u1, u2)
            local_field(1,3) = prefac * local_Ez(a1, a2, b1, b2, u1, u2)
        else
            local_field(1,1) = z_sign * prefac * local_Ex(a1, a2, b1, b2, u1, u2)
            local_field(1,2) = prefac * local_Ey(a1, a2, b1, b2, u1, u2)
            if (dist < 1.0d-12) then
                local_field(1,3) = 1.0d0/(2.0d0*eps0)
            else
                local_field(1,3) = 0.0d0
            end if
        end if

        do i = 1, 3
            field_elem(1,i) = n1(1,i)*local_field(1,1) + N2prime(1,i)*local_field(1,2) + n3(1,i)*local_field(1,3)
        end do

        field_elem = sigma*field_elem
        
    end function field_elem

    function p_calc_noZ(a2, b2, a1, b1, y)
        implicit none
        double precision a1, a2, b1, b2, y, p_calc_noZ, logArg1, logArg2, ans1, ans2
        logArg2 = (1.0d0 + b2*b2)*y + a2*b2 + sqrt((1.0d0 + b2*b2)*y*y + 2*a2*b2*y + a2*a2)
        logArg1 = (1.0d0 + b1*b1)*y + a1*b1 + sqrt((1.0d0 + b1*b1)*y*y + 2*a1*b1*y + a1*a1)
        
        if (logArg2 > 1.0d-14) then
            if (abs(y) > 1.0d-14) then
                ans2 = y * asinh((a2 + b2 * y) / abs(y)) + a2 / sqrt(1.0d0 + b2 * b2) * log(logArg2)
            else
                ans2 = a2 / sqrt(1.0d0 + b2 * b2) * log(logArg2)
            end if
        else
            if (abs(y) > 1.0d-14) then
                ans2 = y * asinh(y * b2 / abs(y))
            else
                ans2 = 0.0d0
            end if
        end if

        if (logArg1 > 1.0d-14) then
            if (abs(y) > 5.0d-14) then
                ans1 = y * asinh((a1 + b1 * y) / abs(y)) + a1 / sqrt(1.0d0 + b1 * b1) * log(logArg1)
            else
                ans1 = a1 / sqrt(1.0d0 + b1 * b1) * log(logArg1)
            end if
        else
            if (abs(y) > 5.0d-14) then
                ans1 = y * asinh(y * b1 / abs(y))
            else
                ans1 = 0.0d0
            end if
        end if

        p_calc_noZ = ans2 - ans1

    end function p_calc_noZ

    function F1(a, b, u)
        implicit none
        double precision F1, a, b, u
        F1 = u*asinh((a + b * u) / sqrt(u * u + 1.0d0))
    end function F1

    function I3(a, b, u1, u2)
        implicit none
        double precision I3, a, b, u1, u2, g1, g2
        g1 = (sqrt(b * b + 1.0d0) * sqrt(a*a + 2.0d0*a*b*u1 + (b*b + 1.0d0)*u1*u1 + 1.0d0) + b*(a + b*u1) + u1)
        g2 = (sqrt(b * b + 1.0d0) * sqrt(a*a + 2.0d0*a*b*u2 + (b*b + 1.0d0)*u2*u2 + 1.0d0) + b*(a + b*u2) + u2)

        if (g1 <= 0.0d0) then
            g1 = -(1.0d0 + a*a + b*b) / (2.0d0*(b*b + 1.0d0)*u1)
        end if
        if (g2 <= 0.0d0) then
            g2 = -(1.0d0 + a*a + b*b) / (2.0d0*(b*b + 1.0d0)*u2)
        end if

        if (abs(g1) < 1.0d-12) then
            if (abs(a) < 1.0d-14) then
                g1 = 1.0d-12
            end if
        end if
        if (abs(g2) < 1.0d-12) then
            if (abs(a) < 1.0d-14) then
                g2 = 1.0d-12
            end if
        end if
                
        I3 = a / sqrt(b*b + 1.0d0) * log(abs(g2 / g1))
    end function I3

    function I4_plus(alpha, gamma, q2, prefac, t1, t2)
        implicit none
        double precision I4_plus, alpha, gamma, q, q2, prefac, t1, t2, g1, g2
        q = sqrt(q2)
        g1 = sqrt(gamma*t1*t1 + alpha)
        g2 = sqrt(gamma*t2*t2 + alpha)

        if (t1 > 1.0d15 .or. t2 > 1.0d15) then
            if (t2 < 1.0d15) then
                I4_plus = (prefac*1.0d0/q*(atan(g2/q) - pi/2.0d0))
                return
            else if (t1 < 1.0d15) then
                I4_plus = (prefac*1.0d0/q*(pi/2.0d0 - atan(g1/q)))
                return
            else
                I4_plus = 0.0d0
                return
            end if
        end if

        I4_plus = prefac*1.0d0/q*atan(q*(g2 - g1)/(q2 + g1*g2))

    end function I4_plus

    recursive function I4(a, b, u1, u2) result(ans)
        implicit none
        double precision a, b, u1, u2, alpha, gamma, q2, prefac, t1, t2, sign, ans
        if (abs(u1 - b/a) < 1.0d-14) then
            if (u2 > b/a) then
                ans = I4(a, b, u1 + 1.0d-14, u2)
                return
            else
                ans = I4(a, b, u1 - 1.0d-14, u2)
                return
            end if
        end if

        if (abs(u2 - b/a) < 1.0d-14) then
            if (u1 > b/a) then
                ans = I4(a, b, u1, u2 + 1.0d-14)
                return
            else
                ans = I4(a, b, u1, u2 - 1.0d-14)
                return
            end if
        end if
        
        if (abs(a) < 1.0d-14) then
            if (a > 0.0d0) then
                ans = I4(a + 1.0d-14, b, u1, u2)
                return
            else
                ans = I4(a - 1.0d-14, b, u1, u2)
                return
            end if
        end if

        alpha = 1.0d0 + (a*a)/(b*b)
        gamma = (a*a + b*b)*(a*a + b*b + 1.0d0)/(b*b)
        q2 = (a*a + b*b)*(a*a + b*b)/(b*b)
        prefac = (a*a/b + b)
        if (a*u1 /= b) then
            t1 = (b*u1 + a)/(a*u1 - b)
        else
            t1 = 1.0d15
        end if
        if (a*u2 /= b) then
            t2 = (b*u2 + a)/(a*u2 - b)
        else
            t2 = 1.0d15
        end if

        sign = 1.0d0

        if (a < 0.0d0) then
            sign = - sign
        end if
        if (b < 0.0d0) then
            sign = -sign
        end if
        if (u1 > b/a) then
            sign = -sign
        end if

        if (((u1<b/a) .and. (b/a<u2)) .or. ((u2<b/a) .and. (b/a<u1))) then
            ans = sign * (I4_plus(alpha, gamma, q2, prefac, t1, 1.0d16) + I4_plus(alpha, gamma, q2, prefac, t2, 1.0d16))
            return
        else
            ans = sign * I4_plus(alpha, gamma, q2, prefac, t1, t2)
            return
        end if
    end function I4

    function I1(a, b, u1, u2)
        implicit none
        double precision I1, a, b, u1, u2
        I1 = F1(a, b, u2) - F1(a, b, u1) + I3(a, b, u1, u2) - I4(a, b, u1, u2)
    end function I1

    function I6(x, u1, u2)
        implicit none
        double precision I6, x, u1, u2
        if (abs(x) < 1.0d-15) then
            I6 = 0.0d0
            return
        end if
        I6 = x*log((sqrt(u2*u2 + x*x + 1.0d0) + u2)/(sqrt(u1*u1 + x*x + 1.0d0) + u1))
    end function I6

    function I7(x, u1, u2)
        implicit none
        double precision I7, x, u1, u2, t1, t2, g1, g2
        if (abs(u1) > 1.0d-16) then
            t1 = 1.0d0/u1
        else
            t1 = 1.0d16
        end if
        if (abs(u2) > 1.0d-16) then
            t2 = 1.0d0/u2
        else
            t2 = 1.0d16
        end if

        g1 = sqrt(1.0d0 + t1*t1*(1.0d0 + x*x))
        g2 = sqrt(1.0d0 + t2*t2*(1.0d0 + x*x))

        I7 = atan(x*(g2 - g1)/(x*x + g2*g1))
    end function I7

    function I2(x, u1, u2)
        implicit none
        double precision I2, x, u1, u2
        if (((u1<0.0d0) .and. (0.0d0<u2)) .or. ((u2<0.0d0) .and. (0.0d0<u1))) then
            if (u1 <= 0.0d0) then
                I2 = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, 0.0d0, abs(u1)) + I7(x, 0.0d0, abs(u2))
            else
                I2 = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, abs(u1), 0.0d0) + I7(x, abs(u2), 0.0d0)
            end if
        else if (u1 <= 0.0d0) then
            I2 = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, u2, u1)
        else
            I2 = (F1(x, 0.0d0, u2) - F1(x, 0.0d0, u1)) + I6(x, u1, u2) + I7(x, u1, u2)
        end if
    end function I2

    function p_calc_o(triangle)         ! 4*pi*eps がどうなっているか確認!!!!!!!!!!!!!!!!!!!!!!
        !
        !
        ! 三角形要素が自身の内心に作る電位の電位係数
        !
        ! Parameters
        ! -----------
        ! triangle : [3,3]の実数行列．それぞれの行に三角形の頂点．
        !            三角形の座標
        !
        ! Returns
        ! -----------
        ! p_calc_o : 実数
        !            電位係数
        !
        !
        implicit none
        double precision triangle(3,3), p_calc_o, a(3), b(3), c(3), ab(3), ac(3), ba(3), bc(3), ca(3), cb(3),&
                        & theta_a, theta_b, theta_c
        a = triangle(1,:); b = triangle(2,:); c = triangle(3,:)
        ab = b-a; ac = c-a; ba = a-b; bc = c-b; ca = a-c; cb = b-c
        theta_a = acos(sum(ab*ac)/(distance(a,b)*distance(a,c)))
        theta_b = acos(sum(bc*ba)/(distance(b,c)*distance(b,a)))
        theta_c = acos(sum(ca*cb)/(distance(c,a)*distance(c,b)))
        p_calc_o = 2.0*(area(triangle)*log10(((1.0+cos(theta_a/2.0))*(1.0+cos(theta_b/2.0))*(1+cos(theta_c/2.0))/&
                & ((1.0-cos(theta_a/2.0))*(1.0-cos(theta_b/2.0))*(1.0-cos(theta_c/2.0))))))/&
                & (distance(a,b)+distance(b,c)+distance(c,a))
    end function p_calc_o

    ! function potential(r, tris, sigmas)
    !     !
    !     !
    !     ! 位置rでの電位を計算する．
    !     ! 
    !     ! Parameters
    !     ! -------------
    !     ! r         : [1,3]の実数配列
    !     !             rの座標
    !     ! tris      : 三角形の構造体を収めた配列
    !     !       `     三角形の座標
    !     ! sigmas    : [n,1]の実数配列．nは三角形の個数．
    !     !             それぞれの三角形の電荷密度
    !     ! 
    !     ! Returns
    !     ! -------------
    !     ! potential : 実数
    !     !             rでの電位
    !     !
    !     !
    !     implicit none
    !     double precision potential, r(1,3), sigmas(:)
    !     integer i ! , n
    !     type(Triangle) tris(:), tri
    !     ! n = ubound(tris, 1)
    !     potential = 0.0d0
    !     do i = 1, num_of_triangles
    !         tri = tris(i)
    !         potential = potential + potential_elem(tri,r,sigmas(i))
    !     end do
    !     potential = potential - sum(E0*r)
    ! end function potential

    function cross_product(a,b)
        implicit none
        double precision a(3), b(3), cross_product(1,3)
        cross_product(1,1) = a(2)*b(3) - a(3)*b(2)
        cross_product(1,2) = a(3)*b(1) - a(1)*b(3)
        cross_product(1,3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product

    function centroid(triangle)
        implicit none
        double precision triangle(3,3), centroid(3)
        centroid = (triangle(1,:) + triangle(2,:) + triangle(3,:))/3
    end function centroid

    function local_Ex(a1, a2, b1, b2, u1, u2)
        implicit none
        double precision local_Ex, a1, a2, b1, b2, u1, u2
        local_Ex = I3p(a2, b2, u1, u2) - I3p(a1, b1, u1, u2)
    end function local_Ex

    function local_Ey(a1, a2, b1, b2, u1, u2)
        implicit none
        double precision local_Ey, a1, a2, b1, b2, u1, u2, I1, I2
        if (abs(b2) > 1.0d-14) then
            I2 = I4_2(a2, b2, u1, u2) + b2*I3p(a2, b2, u1, u2)
        else
            I2 = J2(a2, u1, u2)
        end if

        if (abs(b1) > 1.0d-14) then
            I1 = I4_2(a1, b1, u1, u2) + b1*I3p(a1, b1, u1, u2)
        else
            I1 = J2(a1, u1, u2)
        end if

        local_Ey = I1- I2
        
    end function local_Ey

    function local_Ez(a1, a2, b1, b2, u1, u2)
        implicit none
        double precision local_Ez, a1, a2, b1, b2, u1, u2, I1, I2
        if (abs(b1) > 1.0d-14) then
            I1 = I4(a1, b1, u1, u2);
        else if (((u1<0.0d0) .and. (0.0d0<u2)) .or. ((u2<0.0d0) .and. (0.0d0<u1))) then
            if (u1 <= 0.0d0) then
                I1 = -(I7(a1, 0.0d0, abs(u1)) + I7(a1, 0.0d0, abs(u2)))
            else
                I1 = -(I7(a1, abs(u1), 0.0d0) + I7(a1, abs(u2), 0.0d0))
            end if
        else if (u1 <= 0.0d0) then
            I1 = I7(a1, u1, u2);
        else
            I1 = I7(a1, u2, u1);
        end if

        if (abs(b2) > 1.0d-14) then
            I2 = I4(a2, b2, u1, u2);
        else if (((u1<0.0d0) .and. (0.0d0<u2)) .or. ((u2<0.0d0) .and. (0.0d0<u1))) then
            if (u1 <= 0.0d0) then
                I2 = -(I7(a2, 0.0d0, abs(u1)) + I7(a2, 0.0d0, abs(u2)))
            else
                I2 = -(I7(a2, abs(u1), 0.0d0) + I7(a2, abs(u2), 0.0d0))
            end if
        else if (u1 <= 0.0d0) then
            I2 = I7(a2, u1, u2)
        else
            I2 = I7(a2, u2, u1)
        end if

        local_Ez = I2 - I1;

    end function local_Ez

    function I3p(a, b, u1, u2)
        implicit none
        double precision I3p, a, b, u1, u2, g1, g2
        g1 = (sqrt(b*b + 1.0d0) * sqrt(a*a + 2.0d0*a*b*u1 + (b*b + 1.0d0)*u1*u1 + 1.0d0) + b*(a + b*u1) + u1)
        g2 = (sqrt(b*b + 1.0d0) * sqrt(a*a + 2.0d0*a*b*u2 + (b*b + 1.0d0)*u2*u2 + 1.0d0) + b*(a + b*u2) + u2)
        if (g1 <= 0.0d0) then
            g1 = -(1.0d0 + a*a + b*b)/(2.0d0*(b*b + 1.0d0)*u1)
        end if
        if (g2 <= 0.0d0) then
            g2 = -(1.0d0 + a*a + b*b)/(2.0d0*(b*b + 1.0d0)*u2)
        end if
        I3p = 1.0d0/sqrt(b*b + 1.0d0)*log(abs(g2/g1))
    end function I3p

    function J2(a, u1, u2)
        implicit none
        double precision J2, a, u1, u2, g1, g2
        if (a == 0.0d0) then
            J2 = 0.0d0
        else
            g1 = sqrt(u1*u1 + a*a + 1.0d0)
            g2 = sqrt(u2*u2 + a*a + 1.0d0)
            J2 = a/(2.0d0*abs(a))*log(((g2 - abs(a))*(g1 + abs(a)))/((g2 + abs(a))*(g1 - abs(a))))
        end if
    end function J2

    function I4_2_plus(alpha, gamma, prefac, t1, t2)
        implicit none
        double precision I4_2_plus, alpha, gamma, prefac, t1, t2, g1, g2, q
        g1 = sqrt(alpha * t1 * t1 + gamma)
        g2 = sqrt(alpha * t2 * t2 + gamma)
        if (gamma - alpha <= 0.0d0) then
            I4_2_plus = prefac*(g2 - g1)/((alpha - gamma) + g1*g2)
        else
            q = sqrt(gamma - alpha)
            I4_2_plus = prefac*1.0d0/q*atanh(q*(g2 - g1)/((alpha - gamma) + g1*g2))
        end if
    end function I4_2_plus

    recursive function I4_2(a, b, u1, u2) result(ans)
        implicit none
        double precision  a, b, u1, u2, alpha, gamma, lambda, prefac, t1, t2, ans
        alpha = 1.0d0 + (a*a)/(b*b)
        gamma = (a*a + b*b)*(a*a + b*b + 1.0d0)/(b*b)
        lambda = -a/b
        prefac = (a*a/b + b)
        if (b * u1 /= -a) then 
            t1 = (b - a*u1) / (a + b*u1)
        else
            t1 = 1.0d15
        end if

        if (b * u2 /= -a) then 
            t2 = (b - a*u2) / (a + b*u2)
        else
            t2 = 1.0d15
        end if

        if (((u1<lambda) .and. (lambda<u2)) .or. ((u2<lambda) .and. (lambda<u1))) then
            if (u1 > lambda) then
                ans = I4_2_plus(alpha, gamma, prefac, abs(t1), 1.0d15) + I4_2_plus(alpha, gamma, prefac, abs(t2), 1.0d15)
                return
            else
                ans = I4_2_plus(alpha, gamma, prefac, 1.0d15, abs(t1)) + I4_2_plus(alpha, gamma, prefac, 1.0d15, abs(t2))
                return
            end if
        else if (u1 > lambda) then
            ans = I4_2_plus(alpha, gamma, prefac, t1, t2)
            return
        else
            ans = I4_2_plus(alpha, gamma, prefac, t2, t1)
        end if
    end function I4_2

end subroutine calc_efield