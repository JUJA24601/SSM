import fortran
import calc_matrix
import calc_potential
import numpy as np
from stl import mesh
import time

# 電位計算に解析積分を用いたもの

# 真空誘電率
EPS_0 = 8.8541878128E-12
# 一様電場
E0 = np.array([0.0,0.0,0.0])

def main():
    triangles = mesh.Mesh.from_file("cylinder.stl").vectors     # 三角形のデータ
    
    num_of_triangles = triangles.shape[0]       # 三角形の数
    potential_of_material = 40                 # 物体の表面電位
    
    # 三角形の電荷密度を求めるための準備
    sigmas = np.zeros((num_of_triangles,1))                                 # 電荷密度を収めるための配列
    Vs = np.full((num_of_triangles,1), potential_of_material)               # それぞれの三角形の電位
    P = np.zeros((num_of_triangles,num_of_triangles))                       # 電位係数を収める
    vector = np.zeros((num_of_triangles,1))                                 # 外部電場の影響を考える際に必要になる
    
    # 電位を求めるため，空間を分割
    div = 100       # 空間一辺の分割数
    x = np.linspace(-8,8,div)
    y = np.linspace(-8,8,div)
    z = np.linspace(-8,8,div)
    
    print(f"----- python -----")
    print(f"num_of_triangles : {num_of_triangles}")
    print(f"E0 : {E0}")
    
    start = time.time()         # 開始時間
    
    # ここから物体表面の電荷を求める
    P = calc_matrix.calc_matrix(triangles, num_of_triangles)   # Pを計算
    
    # 外部電位の影響を考慮する (ここもfortranで計算してもいいかも)
    for i in range(num_of_triangles):
        vector[i] = Vs[i] + np.dot(E0, incenter(triangles[i]))
    
    sigmas = 4*np.pi*EPS_0*np.linalg.solve(P, vector)       # 電荷密度を求める (ここもfortranでしてもいいかも?)
    
    mid = time.time()
    time_diff1 = mid - start
    print(f"time 1 : {time_diff1}")
    
    # 電荷密度を保存
    np.save("sigmas", sigmas)   # 保存するデータの名前を変更
    # np.save("P_numerical", P/(4*np.pi*eps_0))
    
    # ここから空間の電位を求める
    Potential = calc_potential.calc_potential(triangles, sigmas, [E0], x, y, z, num_of_triangles, div)          # 空間内の電位を求める
    np.save("Potential", Potential)     # 空間内の電位を保存
    
    end = time.time()       # 終了時間
    
    time_diff1 = mid - start
    print(f"time 1 : {time_diff1}")
    time_diff2 = end - mid
    print(f"time 2 : {time_diff2}")
    time_diff = end - start
    print(f"elapsed time : {time_diff}")

def incenter(triangle) :
    """
    三角形の内心を求める
    
    Parameters
    ----------
    triangle : float (3,3)
        三角形の座標, triangle.shape -> (3,3) の必要あり
    
    Returns
    ---------
    incentre : float (1,3)
        三角形の内心
    """
    a = triangle[0]
    b = triangle[1]
    c = triangle[2]
    ab = distance(a,b)
    bc = distance(b,c)
    ca = distance(c,a)
    incentre = (bc*a+ca*b+ab*c)/(ab+bc+ca)
    return incentre

def distance(a, b) :
    """
    2点間 a, b の距離を求める
    
    Parameters
    ----------
    a : float (3,1)
        aの座標
    b : float (3,1)
        bの座標
    
    Returns
    ----------
    distance : float
        2点間の距離
    """
    distance = np.linalg.norm(b-a)
    
    return distance

if __name__ == '__main__':
    main()