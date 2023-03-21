# -*- coding: utf-8 -*-#
# lnlamp内部使用的python脚本,其核心是进行多项式拟合
import sys
import matplotlib.pyplot as plt
import numpy as np


def matrix_n(alist=[],n=0):
    llen = len(alist)
    r_list = [0]*llen
    for i in range(llen):
        r_list[i] = alist[i] ** n
    return r_list 


# 将结果保存为txt文件
def writer_txt(file, nd_name):
    np.savetxt(file, nd_name, delimiter="\t", fmt="%s")



# # 检查python文件参数.
# if len(sys.argv) > 19:
#     print('参数个数为:', len(sys.argv), '个参数。')
#     print('参数列表:', str(sys.argv))
#     print('脚本名为：', sys.argv[0])
#     for i in range(1, len(sys.argv)):
#         print('参数 %s 为：%s' % (i, sys.argv[i]))
#     print("lnlamp: error: lnlamp_internal_curve_fitting.py 接收的参数存在错误!!!")
#     exit()

# 解析输入
x=[]
y=[]
for i in range(1, len(sys.argv)):
    if i%2 == 1 :
        x.append(float(sys.argv[i]))
    else:
        y.append(float(sys.argv[i]))
# 检查获取的点是一一对应的
print("rate(%):")
print(x)
print("time(s):")
print(y)
if len(x) != len(y):
    print("lnlamp: error: lnlamp_internal_curve_fitting.py 接收的参数不是一一对应的，存在问题!!!")
    exit() 
m = len(x)

# 初始化矩阵
x_0 = [1]*m
x_2 = matrix_n(x,2)
x_3 = matrix_n(x,3)

x_matrix = np.matrix((x_0,x,x_2),dtype=np.float64)
# print(x_matrix)

y_matrix = np.matrix((y),dtype=np.float64)
y_matrix = y_matrix.T

x_matrix_t = x_matrix.T
# print(x_matrix_t)

x_end=x_matrix * x_matrix_t
b_end=x_matrix * y_matrix
print("X系数矩阵:")
print(x_end)
print("B常数列矩阵:")
print(b_end)
out = ((x_matrix * x_matrix_t)** -1)*(x_matrix * y_matrix)
extreme_value = (-1.0 * float(out[1]) ) / (2.0 * float(out[2]))
print("解向量A,也即(AX=B)的解向量:")
print(out)

print("二次拟合曲线的极值点:")
print(extreme_value)

# 输出原始的点
plt.plot(x,y,"ko")
# 输出拟合的曲线
xx = np.linspace(x[0], x[len(x)-1])
yy = float(out[2])*xx*xx + float(out[1])*xx + float(out[0])
plt.plot(xx,yy,label="time="+"{:.{}e}".format(float(out[2]), 3)+"*r^2+ "+"{:.{}e}".format(float(out[1]), 3)+"*r+ "+"{:.{}e}".format(float(out[0]), 5),color="red")

# 设置横纵坐标轴和图例
plt.ylabel('time(s)', fontsize=10.5)
plt.xlabel('rate', fontsize=10.5)
plt.legend()
# 控制纵坐标显示科学计算法，并显示三个有效数字
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText=True)
plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.2e'))

# 调整图像布局
plt.tight_layout()

# 保存结果 or 暂时结果
# plt.show()
plt.savefig("lnlamp_predict_result.png")
out = np.insert(out, 3, extreme_value, 0) # 第0维度(行)第3行添加[1,1,1,1]
writer_txt("lnlamp_temp_result.txt",out)
exit(0)
