import numpy as np
import sympy
import matplotlib
import matplotlib.pyplot as plt

# solve b =Ax
# np.set_printoptions(np.inf)

def solve():
    path = "data0805/"

    file_Areal = path + "data_A_real.csv"
    file_Aimg = path + "data_A_img.csv"
    A_real = np.loadtxt(file_Areal,delimiter=",")
    A_img = np.loadtxt(file_Aimg,delimiter=",")
    A = A_real + A_img * 1.0j

    file_Umea_real = path + "data_Umea_real.csv"
    file_Umea_img = path + "data_Umea_img.csv"
    Umea_real = np.loadtxt(file_Umea_real,delimiter=",")
    Umea_img = np.loadtxt(file_Umea_img,delimiter=",")
    Umea = Umea_real + Umea_img * 1.0j
    print(Umea)
    # temp = np.abs(Umea)
    # print_matrix(temp)
    Umea_db = set_magdb(Umea)
    # print_matrix(result)

    file_ans_real = path + "data_ans_real.csv"
    file_ans_img = path + "data_ans_img.csv"
    ans_real = np.loadtxt(file_ans_real,delimiter=",")
    ans_img = np.loadtxt(file_ans_img,delimiter=",")
    ans = ans_real + ans_img * 1.0j
    result_db = set_magdb(np.dot(A,ans))

    print(result_db.shape)
    print(Umea_db.shape)
    # y = np.stack([Umea_db,result_db])
    y = [Umea_db,result_db]
    print(y)
    print_matrix(y)


    # b = np.array([3+4j,6+0j,8+6j])
    # ans = np.linalg.solve(A, b)
    # print("ans x =",ans)
    # print("residual error =",b - np.dot(A,ans))
    # print("ans Ax =",np.dot(A,ans))

# グラフ表示
def print_matrix(y):
    row = len(y[0])
    print("row=",row)

    x = np.linspace(1,row+1,row)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y[0])

    for i in range(len(y)):
        ax.plot(x,y[i],"label"+str(i))

    plt.legend()
    ax.grid()
    plt.show()

# Umea -> U[db]
def set_magdb(Umea,val_max = 0):
    assert Umea.shape[0]%2 == 0 , "Umea size is wrong"
    row = int(Umea.shape[0] / 2)
    data = np.sqrt(np.square(np.abs(Umea[0:row])) + np.square(np.abs(Umea[row:2*row])))
    val_max = np.amax(data)
    print(val_max)
    result = 20 * np.log10(data/val_max)

    return result


solve()
