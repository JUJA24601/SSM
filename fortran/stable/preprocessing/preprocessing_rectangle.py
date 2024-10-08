import numpy as np

def main():
    # ここで四角形を作成する
    div = 100
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
            
    num_of_rectangles = rectangles.shape[0]        # 三角形の数
    
    header = "apex1_x, apex1_y, apex1_z, apex2_x, apex2_y, apex2_z, apex3_x, apex3_y, apex3_z, apex4_x, apex4_y, apex4_z, a, b, p0_x, p0_y, p0_z, n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, n3_x, n3_y, n3_z"
    csv = np.empty((num_of_rectangles, 26))
    
    for i in range(num_of_rectangles) :
        coord = rectangles[i]
        a = distance(coord[0], coord[1])
        b = distance(coord[0], coord[3])
        p0 = coord[0]
        n1 = (coord[1] - coord[0])/a
        n2 = (coord[3] - coord[0])/b
        n3 = np.cross(n1,n2)
        temp = np.concatenate([coord.ravel(), [a], [b], p0, n1, n2, n3])
        csv[i] = np.copy(temp)
        
    np.savetxt('data.csv', csv, delimiter=',', header=header, comments='')

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