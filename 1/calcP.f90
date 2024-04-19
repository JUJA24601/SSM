subroutine calc_p(triangles, N, num, P)
    implicit none
    integer i, j, N, num
    double precision triangles(num,3,3), P(num,num), eps0, pi, E0(3)
    intent(in) triangles, N, num
    intent(out) P

    eps0 = 8.8541878128d-12     ! 誘電率
    pi = acos(-1.0d0)           ! 円周率
    ! E0 = [20.0d0,0.0d0,0.0d0]   ! 背景電界
    print *, "----- fortran (calcP) -----"
    print *, "num : ",num
    print *, "N : ", N
    ! print *, "E0 : ", E0

    do i = 1,num
        do j = 1,num
            if (i /= j) then
                P(j,i) = p_calc(triangles(i,:,:), triangles(j,:,:), N)
            else
                P(j,i) = p_calc_o(triangles(i,:,:))
            end if
        end do
    end do

contains

    function potential_elem(triangle, r, N, sigma)
        !
        !
        ! 三角形要素triangleによってできる任意の点rでの電位を求める．
        ! 
        ! Parameters
        ! ---------
        ! triangle : [3,3]の実数配列でそれぞれの行に頂点が保存される．
        !            三角形の座標
        ! r        : [1,3]の実数配列
        !            任意の点
        ! N        : 整数               <-------------------- 削除の可能性
        !            三角形の一辺の分割数
        ! sigma    : 実数
        !            三角形の電荷密度
        !
        ! Returns
        ! ---------
        ! potential_elem : 実数
        !                  三角形要素によって点rで生じる電位
        ! 
        !
        implicit none
        double precision triangle(3,3), r(3), sigma, potential_elem, &
                        & d, w, point(3), a(3), b(3), c(3), ab(3), bc(3)
        integer N, i, j
        potential_elem = 0.0
        a = triangle(1,:); b = triangle(2,:); c = triangle(3,:)
        ab = b-a; bc = c-b
        d = 1.0/(3*N*N)
        do i = 0, N
            do j = 0, i
                point = a + i*ab/N + j*bc/N
                if ((i==0 .or. i==N) .and. (j==0 .or. j==N)) then
                    w = d
                else if ((i>0 .and. i<N .and. (j==0 .or. j==i)) .or. (i==N .and. j>0 .and. j<N)) then
                    w = 3*d
                else
                    w = 6*d
                endif
                potential_elem = potential_elem + w/distance(point, r)
            enddo
        enddo
        potential_elem = (potential_elem*sigma*area(triangle))/(4*pi*eps0)
    end function potential_elem


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


    function p_calc(triangle_j, triangle_i, N)
        ! 
        !
        ! 電位係数を求める．(三角形jが三角形iの内心に作る)
        !
        ! Parameters
        ! -----------
        ! triangle_j : [3,3]の実数配列でそれぞれの行に三角形の頂点が保存される
        !              三角形jの座標
        ! triangle_i : [3,3]の実数配列でそれぞれの行に三角形の頂点が保存される
        !              三角形iの座標
        ! N          : 整数                     <-------------------------削除するかも
        !              三角形の一辺の分割数
        !
        ! Returns
        ! -----------
        ! p_calc : 実数
        !          電位係数
        !
        !
        implicit none
        double precision triangle_i(3,3), triangle_j(3,3), p_calc, d, w, point(3), a(3), b(3), c(3), ab(3), bc(3)
        integer N, i, j
        p_calc = 0.0
        a = triangle_j(1,:); b = triangle_j(2,:); c = triangle_j(3,:)
        ab = b-a; bc = c-b
        d = 1.0/(3*N*N)
        do i = 0, N
            do j = 0, i
                point = a + i*ab/N + j*bc/N
                if ((i==0 .or. i==N) .and. (j==0 .or. j==N)) then
                    w = d
                else if ((i>0 .and. i<N .and. (j==0 .or. j==i)) .or. (i==N .and. j>0 .and. j<N)) then
                    w = 3.0*d
                else
                    w = 6.0*d
                endif
                p_calc = p_calc + w/distance(point, incenter(triangle_i))
            enddo
        enddo
        p_calc = (p_calc*area(triangle_j))
    end function p_calc

    function p_calc_o(triangle)
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

    function potential(r, triangles, sigmas)
        !
        !
        ! 位置rでの電位を計算する．
        ! 
        ! Parameters
        ! -------------
        ! r         : [1,3]の実数配列
        !             rの座標
        ! triangles : [n,3,3]の実数配列．nは三角形の個数．
        !       `     三角形の座標
        ! sigmas    : [n,1]の実数配列．nは三角形の個数．
        !             それぞれの三角形の電荷密度
        ! 
        ! Returns
        ! -------------
        ! potential : 実数
        !             rでの電位
        !
        !
        implicit none
        double precision potential, triangles(:,:,:), r(3), sigmas(:), triangle(3,3)
        integer i, n
        n = ubound(triangles, 1)
        potential = 0.0d0
        do i = 1, n
            triangle = triangles(i,:,:)
            potential = potential + potential_elem(triangle,r,N,sigmas(i))
        end do
        potential = potential - sum(E0*r)
    end function potential

end subroutine calc_p