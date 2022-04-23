import math
import numpy as np

def coordinate(pathway):
    #读取文件
    f = open(pathway)
    poscar = f.read()
    f.close()
    #数据除噪
    poscar = ' '.join(poscar.split())
    poscar = poscar.replace(" T T T ","").replace(" F F F ","").replace(" ", ",")
    poscar = poscar.split(",")

    for i in range(0,14):
        del poscar[0]
        i =+ 1

    for i in range(len(poscar)):
        poscar[i] = float(poscar[i])

    #转换成坐标
    atom = ["c"+str(i) for i in range(len(poscar)//3)]

    for i in range(len(poscar)//3):
        atom[i] = np.array([poscar[3*i],poscar[3*i+1],poscar[3*i+2]])

    return atom

def poav(A, B, C, O):
    #计算sigma键矢量
    OA = O - A
    OB = O - B
    OC = O - C

    #计算sigma键模长
    OA_l = np.sqrt(OA.dot(OA))
    OB_l = np.sqrt(OB.dot(OB))
    OC_l = np.sqrt(OC.dot(OC))

    #计算POAV步骤
    a1 = OB_l*(A[0] - O[0]) - OA_l*(B[0] - O[0])
    a2 = OB_l*(A[1] - O[1]) - OA_l*(B[1] - O[1])
    a3 = OB_l*(A[2] - O[2]) - OA_l*(B[2] - O[2])

    b1 = OC_l*(A[0] - O[0]) - OA_l*(C[0] - O[0])
    b2 = OC_l*(A[1] - O[1]) - OA_l*(C[1] - O[1])
    b3 = OC_l*(A[2] - O[2]) - OA_l*(C[2] - O[2])

    c1 = OC_l*(B[0] - O[0]) - OB_l*(C[0] - O[0])
    c2 = OC_l*(B[1] - O[1]) - OB_l*(C[1] - O[1])
    c3 = OC_l*(B[2] - O[2]) - OB_l*(C[2] - O[2])

    l = (a2*b3-a3*b2)/math.sqrt((a2*b3-a3*b2)**2 + (a3*b1-a1*b3)**2 + (a1*b2-a2*b1)**2)
    m = (a3*b1-a1*b3)/math.sqrt((a2*b3-a3*b2)**2 + (a3*b1-a1*b3)**2 + (a1*b2-a2*b1)**2)
    n = (a1*b2-a2*b1)/math.sqrt((a2*b3-a3*b2)**2 + (a3*b1-a1*b3)**2 + (a1*b2-a2*b1)**2)

    poa_r = np.arccos(abs(l*(A[0]-O[0]) + m*(A[1]-O[1]) + n*(A[2]-O[2]))/math.sqrt((A[0]-O[0])**2 + (A[1]-O[1])**2 + (A[2]-O[2])**2))
    poa = poa_r * 180 / np.pi
    poav = 90 - poa

    return poav

#在此输入POSCAR文件路径
pathway = "C:/Users/Jarvis/Desktop/POSCAR.vasp"
atom = coordinate(pathway)

atoms = ["A", "B", "C", "O"]

for j in range(len(atom)):
    num = 0
    for i in range(len(atom)):
        #通过原子间距筛选ABCO四个原子
        bond = atom[j] - atom[i]
        bond_l = np.sqrt(bond.dot(bond))
        #筛选中心原子并记录编号
        if bond_l ==0:
            atoms[3] = atom[i]
            site_0 = atom[i]
        #筛选ABC原子并记录编号
        elif  0 < bond_l < 1.5:
            atoms[num] = atom[i]
            num += 1
            if num == 1:
                site_1 = i + 1
            if num == 2:
                site_2 = i + 1
            if num == 3:
                site_3 = i + 1
    #剔除边缘原子（无法求POAV）
    if num == 3:
        #打印中心原子及周围三个原子编号、中心原子坐标、中心原子的POAV
        print(j + 1, site_1, site_2, site_3, site_0, poav(atoms[0], atoms[1], atoms[2], atoms[3]))


bonds = []
#计算任意两原子间距
for j in range(len(atom)):
    for i in range(len(atom)):
        bond = atom[j] - atom[i]
        bond_l = np.sqrt(bond.dot(bond))
        #筛选相邻两原子键长
        if 0 < bond_l < 1.5:
            bonds.append(bond_l)

#去重
new_bonds = []
for i in bonds:
    if i not in new_bonds:
        new_bonds.append(i)

new_bonds.sort()
sum = 0
#求和
for i in new_bonds:
    sum += i
#计算平均键长
avg = sum/len(new_bonds)
print("平均键长：", avg)
