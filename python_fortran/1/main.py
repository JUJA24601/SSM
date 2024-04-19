import calcP
import numpy as np
from stl import mesh
import time

# 数値積分を用いる

# 真空誘電率
eps_0 = 8.8541878128E-12
# 三角形の一辺の分割数
N = 1
# 一様電場
E0 = np.array([0.0,0.0,0.0])

def main():
    triangles = mesh.Mesh.from_file("sphere.stl").vectors     # 取り込む三角形データの名前に変更
    
    num = triangles.shape[0]        # 三角形の数
    sigmas = np.zeros((num,1))      # 電荷密度を収めるための配列
    Vs = np.full((num,1), 40)       # それぞれの三角形の電位
    P = np.zeros((num,num))         # 電位係数を収める
    vector = np.zeros((num,1))
    
    # 電位を求めるため，空間を分割
    div_x = 50       # 空間一辺の分割数
    div_y = 50       # 空間一辺の分割数
    div_z = 50       # 空間一辺の分割数
    x = np.linspace(-8,8,div_x)
    y = np.linspace(-8,8,div_y)
    z = np.linspace(-8,8,div_z)
    
    print(f"----- python -----")
    print(f"num : {num}")
    print(f"N : {N}")
    print(f"E0 : {E0}")
    
    start = time.time()         # 開始時間


    # ここから物体表面の電荷を求める
    P = calcP.calc_p(triangles, N, num)   # Pを計算

    # ここで電場による影響を考慮する
    for i in range(num):
        vector[i] = Vs[i] + np.dot(E0, incenter(triangles[i]))

    sigmas = 4*np.pi*eps_0*np.linalg.solve(P, vector)

    mid = time.time()
    time_diff1 = mid - start
    print(f"time 1 : {time_diff1}")

    np.save("sigmas", sigmas)
    
    
    # ここから空間の電位を求める
    Potential = calcV.calc_v(triangles, sigmas, E0, x, y, z, N, num, div_x, div_y, div_z)          # 空間内の電位を求める
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
    triangle : ndarray
        三角形の座標, triangle.shape -> (3,3) の必要あり

    Returns
    ---------
    incentre : ndarray
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
    a : darray
        aの座標
    b : darray
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