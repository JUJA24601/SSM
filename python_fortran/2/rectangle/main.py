import calc_matrix
import numpy as np
from stl import mesh
import time

# 電位計算に解析積分を用いたもの

# 真空誘電率
EPS_0 = 8.8541878128E-12
# 一様電場
E0 = np.array([0.0,0.0,0.0])

def main():
    # triangles = mesh.Mesh.from_file("cylinder.stl").vectors     # 三角形のデータ
    
    # # ここで立方体を作る
    # div = 6
    # for i in range(div):
    #     for j in range(div):
    #         temp = np.array([[,,],
    #                          [,,],
    #                          [,,],
    #                          [,,]])
            
    # ここで平面四角形を作成する
    div = 49
    rectangles = np.zeros(((2*div)**2,4,3))
    print(rectangles.shape)
    tik = np.zeros(2*div+1)
    tik[0] = -0.5
    tik[2*div] = 0.5
    temp = 0
    for i in range(div):
        r_0 = 0.5*(1 - (i/div)**2)
        r_1 = 0.5*(1 - ((i+1)/div)**2)
        temp += (r_0 - r_1)
        tik[2*div-1-i] = 0.5 - temp
        tik[i+1] = -(0.5 - temp)
        # print(0.5 - temp)
    # print(temp)
    # print(tik)
    count = 0
    for i in range(2*div):
        x_1 = tik[i]
        x_2 = tik[i+1]
        # print(x_1, x_2)
        for j in range(2*div):
            y_1 = tik[j]
            y_2 = tik[j+1]
            temp = np.array([[x_1,y_1,0],
                             [x_2,y_1,0],
                             [x_2,y_2,0],
                             [x_1,y_2,0]])
            rectangles[count] = temp
            count += 1
    
    # print(rectangles[0])
    # print(rectangles[1])
    # print(rectangles[(2*div)**2-1])
    
            
    # rectangles = np.array([
    #     [[-0.5,-0.5, 0],
    #      [ 0.5,-0.5, 0],
    #      [ 0.5, 0.5, 0],
    #      [-0.5, 0.5, 0]]
    # ])
    
    num_of_rectangles = rectangles.shape[0]       # 4角形の数
    potential_of_material = 40                 # 物体の表面電位
    
    # 4角形の電荷密度を求めるための準備
    sigmas = np.zeros((num_of_rectangles,1))                                 # 電荷密度を収めるための配列
    Vs = np.full((num_of_rectangles,1), potential_of_material)               # それぞれの4角形の電位
    P = np.zeros((num_of_rectangles,num_of_rectangles))                       # 電位係数を収める
    # vector = np.zeros((num_of_triangles,1))                                 # 外部電場の影響を考える際に必要になる
    
    print(f"----- python -----")
    print(f"num_of_rectangles : {num_of_rectangles}")
    # print(f"E0 : {E0}")
    
    start = time.time()         # 開始時間
    
    # ここから物体表面の電荷を求める
    P = calc_matrix.calc_matrix(rectangles, num_of_rectangles)   # Pを計算
    print(P.shape)
    
    # # 外部電場の影響を考慮する (ここもfortranで計算してもいいかも)
    # for i in range(num_of_triangles):
    #     vector[i] = Vs[i] + np.dot(E0, incenter(triangles[i]))
    
    sigmas = 4*np.pi*EPS_0*np.linalg.solve(P, Vs)       # 電荷密度を求める (ここもfortranでしてもいいかも?)
    # for i in range(num_of_rectangles):
    #     print(sigmas[i])
    
    mid = time.time()
    time_diff1 = mid - start
    # print(f"time 1 : {time_diff1}")
    
    # # 電荷密度を保存
    np.save("sigmas", sigmas)   # 保存するデータの名前を変更
    # np.save("P_numerical", P/(4*np.pi*eps_0))
    
    # ここから静電容量を求める
    Q = 0
    for i in range(num_of_rectangles) : 
        Q += sigmas[i]*area(rectangles[i])
    # print(Q)
    
    C = Q/potential_of_material
    print(f"Capacitance : {C}")
    print(f"Capacitance/4pi*eps0 : {C/(4*np.pi*EPS_0)}")
    
    # # ここから空間の電位を求める
    # Potential = calc_potential.calc_potential(triangles, sigmas, [E0], x, y, z, num_of_triangles, div)          # 空間内の電位を求める
    # np.save("Potential", Potential)     # 空間内の電位を保存
    
    end = time.time()       # 終了時間
    
    time_diff1 = mid - start
    # print(f"time 1 : {time_diff1}")
    time_diff2 = end - mid
    # print(f"time 2 : {time_diff2}")
    time_diff = end - start
    # print(f"elapsed time : {time_diff}")

def area(rectangle) :
    a = distance(rectangle[0], rectangle[1])
    b = distance(rectangle[0], rectangle[3])
    area = a*b
    return area

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