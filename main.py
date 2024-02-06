import math
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from openpyxl.workbook import Workbook

def funp3(x):
    b = x ** 4 - 8 * x ** 3 + 86 * x ** 2 - 280 * x + 80
    d = (x ** 8 - 16 * x ** 7 + 196 * x ** 6 - 1456 * x ** 5 + 9796 * x **4 - 44320 * x ** 3 + 87040 * x ** 2 - 44800 * x + 6400)
    a = 40 * (x ** 2 - 4 * x + 8)
    return a, b, d

def koeff(x1, p3):
    p1 = 2
    p2 = 6
    p4 = 10 * p3
    x2 = (p2 * x1 + (x1 - p1) * (1 + 2 * p3)) / x1 ** 2
    x4 = x2 + ((x1 - p1) * (1 + 2 * p3)) / p4
    x3 = 2 * p1 - x1
    return x2, x3, x4

def jacobi(x1, x2, x3, x4, p3):
    p2 = 6
    p4 = 10 * p3
    j = np.zeros((4, 4))
    j[0][0] = -p2 - 1 + 2 * x1 * x2 - p3
    j[0][1] = x1 ** 2
    j[0][2] = p3
    j[0][3] = 0
    j[1][0] = p2 - 2 * x1 * x2
    j[1][1] = -(x1 ** 2) - p4
    j[1][2] = 0
    j[1][3] = p4
    j[2][0] = p3
    j[2][1] = 0
    j[2][2] = -p2 - 1 + 2 * x3 * x4 - p3
    j[2][3] = x3 ** 2
    j[3][0] = 0
    j[3][1] = p4
    j[3][2] = p2 - 2 * x3 * x4
    j[3][3] = -(x3 ** 2) - p4
    eigvals = np.linalg.eigvals(j)
    return eigvals, j

def examfun(x1, x2, x3, x4, p3):
    dx1 = 2 - 7 * x1 + x1 ** 2 * x2 + p3 * (x3 - x1)
    dx2 = 6 * x1 - x1 ** 2 * x2 + 10 * p3 * (x4 - x2)
    dx3 = 2 - 7 * x3 + x3 ** 2 * x4 - p3 * (x3 - x1)
    dx4 = 6 * x3 - x3 ** 2 * x4 - 10 * p3 * (x4 - x2)
    return dx1, dx2, dx3, dx4

def solvex():
    np.set_printoptions(precision=3)
    h = 0.05
    x0 = 0.2
    x_end = 4
    n = int((x_end - x0) / h)
    x1 = np.zeros(n)
    x2_1 = np.zeros(n)
    x3_1 = np.zeros(n)
    x4_1 = np.zeros(n)
    x2_2 = np.zeros(n)
    x3_2 = np.zeros(n)
    x4_2 = np.zeros(n)
    j1 = np.zeros((n, 4), dtype=np.complex_)
    j1_sol = np.zeros(n)
    j2_sol = np.zeros(n)
    j2 = np.zeros((n, 4), dtype=np.complex_)
    p31 = np.zeros(n)
    p32 = np.zeros(n)
    b1 = []
    b2 = []
    jacobi1 = np.zeros((n, 4, 4))
    jacobi2 = np.zeros((n, 4, 4))
    for i in range(n):
        x1[i] = x0 + h * i
    for i in range(n):
        j1_sol[i] = -1
        j2_sol[i] = -1
        a, b, d = funp3(x1[i])
        if d >= 0 and a != 0:
            p31[i] = (-b + math.sqrt(d)) / a
            x2_1[i], x3_1[i], x4_1[i] = koeff(x1[i], p31[i])
            j1[i], jacobi1[i] = jacobi(x1[i], x2_1[i], x3_1[i], x4_1[i], p31[i])
            j1_sol[i] = max(j1[i]) < 0
            if j1_sol[i] != j1_sol[i - 1] and i != 0 and j1_sol[i - 1] != -1 and j1_sol[i] != -1:
                b1.append(x1[i])
            p32[i] = (-b - math.sqrt(d)) / a
            x2_2[i], x3_2[i], x4_2[i] = koeff(x1[i], p32[i])
            j2[i], jacobi2[i] = jacobi(x1[i], x2_2[i], x3_2[i], x4_2[i], p32[i])
            j2_sol[i] = max(j2[i]) < 0
            if j2_sol[i] != j2_sol[i - 1] and i != 0 and j2_sol[i - 1] != -1 and j2_sol[i] != -1:
                b2.append(x1[i])
        else:
            p31[i] = 0
            p32[i] = 0
    return x1, p31, p32, j1, j2, x2_1, x2_2, x3_1, x3_2, x4_1, x4_2, j1_sol, j2_sol, b1, b2, n, jacobi1, jacobi2

def examblock(x1, x2_1, x2_2, x3_1, x3_2, x4_1, x4_2, p31, p32, n):
    dx1_1 = np.zeros(n)
    dx2_1 = np.zeros(n)
    dx3_1 = np.zeros(n)
    dx4_1 = np.zeros(n)
    dx1_2 = np.zeros(n)
    dx2_2 = np.zeros(n)
    dx3_2 = np.zeros(n)
    dx4_2 = np.zeros(n)
    for i in range(n):
        dx1_1[i], dx2_1[i], dx3_1[i], dx4_1[i] = examfun(x1[i], x2_1[i], x3_1[i], x4_1[i], p31[i])
        dx1_2[i], dx2_2[i], dx3_2[i], dx4_2[i] = examfun(x1[i], x2_2[i], x3_2[i], x4_2[i], p32[i])
    return dx1_1, dx1_2, dx2_1, dx2_2, dx3_1, dx3_2, dx4_1, dx4_2

print("---------------------MAIN---------------------------")

x1, p31, p32, j1, j2, x2_1, x2_2, x3_1, x3_2, x4_1, x4_2, j1_sol, j2_sol, b1, b2, n, jacobi1, jacobi2 = solvex()
print(b1, b2)
zip(x1, p31, x2_1, x3_1, x4_1)
with open('first_solution.csv', 'w') as f:
    writer = csv.writer(f)
    header = ['x1', 'p3', 'x2', 'x3', 'x4']
    writer.writerow(header)
    writer.writerows(zip(x1, p31, x2_1, x3_1, x4_1))
zip(x1, p32, x2_2, x3_2, x4_2)
with open('second_solution.csv', 'w') as f:
    writer = csv.writer(f)
    header = ['x1', 'p3', 'x2', 'x3', 'x4']
    writer.writerow(header)
    writer.writerows(zip(x1, p32, x2_2, x3_2, x4_2))
zip(x1, p31, p32, j1, j2, j1_sol, j2_sol)
with open('jacobi.csv', 'w') as f:
    writer = csv.writer(f)
    header = ['x1', 'p31', 'p32', 'j1', 'j2', '1 sol', '2 sol']
    writer.writerow(header)
    writer.writerows(zip(x1, p31, p32, j1, j2, j1_sol, j2_sol))
first_solution = pd.read_csv('first_solution.csv')
second_solution = pd.read_csv('second_solution.csv')
jacobi_solution = pd.read_csv('jacobi.csv')
with pd.ExcelWriter('solution.xlsx') as writer:
    first_solution.to_excel(writer, sheet_name='first_solution')
    second_solution.to_excel(writer, sheet_name='second_solution')
    jacobi_solution.to_excel(writer, sheet_name='jacobi_solution')
print("-------------------Examination------------------")
dx1_1, dx1_2, dx2_1, dx2_2, dx3_1, dx3_2, dx4_1, dx4_2 = examblock(x1, x2_1,
x2_2, x3_1, x3_2, x4_1, x4_2, p31, p32, n)
zip(dx1_1, dx2_1, dx3_1, dx4_1)
with open('exam1.csv', 'w') as f:
    writer = csv.writer(f)
    header = ['dx1', 'dx2', 'dx3', 'dx4']
    writer.writerow(header)
    writer.writerows(zip(dx1_1, dx2_1, dx3_1, dx4_1))
zip(dx1_2, dx2_2, dx3_2, dx4_2)
with open('exam2.csv', 'w') as f:
    writer = csv.writer(f)
    header = ['dx1', 'dx2', 'dx3', 'dx4']
    writer.writerow(header)
    writer.writerows(zip(dx1_2, dx2_2, dx3_2, dx4_2))
exam1 = pd.read_csv('exam1.csv')
exam2 = pd.read_csv('exam2.csv')
with pd.ExcelWriter('examination.xlsx') as writer:
    exam1.to_excel(writer, sheet_name='exam1')
    exam2.to_excel(writer, sheet_name='exam2')
for i in range(n):
    if j1[i].max() < 0 and p31[i] > 0.0:
        plt.scatter(p31[i], x1[i], c='red', s=1)
    elif p31[i] > 0.0:
        plt.scatter(p31[i], x1[i], c='blue', s=1)
for i in range(n):
    if j2[i].max() < 0 and p32[i] > 0.0:
        plt.scatter(p32[i], x1[i], c='red', s=1)
    elif p32[i] > 0.0:
        plt.scatter(p32[i], x1[i], c='blue', s=1)
plt.title("Зависимость x1 от p3")
plt.xlabel("p3")
plt.ylabel("x1")
plt.show()
for i in range(n):
    if j1[i].max() < 0 and p31[i] > 0.0:
        plt.scatter(p31[i], x2_1[i], c='red', s=1)
    elif p31[i] > 0.0:
        plt.scatter(p31[i], x2_1[i], c='blue', s=1)
for i in range(n):
    if j2[i].max() < 0 and p32[i] > 0.0:
        plt.scatter(p32[i], x2_2[i], c='red', s=1)
    elif p32[i] > 0.0:
        plt.scatter(p32[i], x2_2[i], c='blue', s=1)
plt.title("Зависимость x2 от p3")
plt.xlabel("p3")
plt.ylabel("x2")
plt.show()
for i in range(n):
    if j1[i].max() < 0 and p31[i] > 0.0:
        plt.scatter(p31[i], x3_1[i], c='red', s=1)
    elif p31[i] > 0.0:
        plt.scatter(p31[i], x3_1[i], c='blue', s=1)
for i in range(n):
    if j2[i].max() < 0 and p32[i] > 0.0:
        plt.scatter(p32[i], x3_2[i], c='red', s=1)
    elif p32[i] > 0.0:
        plt.scatter(p32[i], x3_2[i], c='blue', s=1)
plt.title("Зависимость x3 от p3")
plt.xlabel("p3")
plt.ylabel("x3")
plt.show()
for i in range(n):
    if j1[i].max() < 0 and p31[i] > 0.0:
        plt.scatter(p31[i], x4_1[i], c='red', s=1)
    elif p31[i] > 0.0:
        plt.scatter(p31[i], x4_1[i], c='blue', s=1)
for i in range(n):
    if j2[i].max() < 0 and p32[i] > 0.0:
        plt.scatter(p32[i], x4_2[i], c='red', s=1)
    elif p32[i] > 0.0:
        plt.scatter(p32[i], x4_2[i], c='blue', s=1)
plt.title("Зависимость x4 от p3")
plt.xlabel("p3")
plt.ylabel("x4")
plt.show()
