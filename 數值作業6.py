# 題目一：Gaussian elimination with pivoting
A1 = [
    [1.19, 2.11, -100, 1],
    [14.2, -0.112, 12.2, -1],
    [0, 100, -99.9, 1],
    [15.3, 0.110, -13.1, -1]
]
b1 = [1.12, 3.44, 2.15, 4.16]

# 高斯消去法
n = 4
for i in range(n):
    max_row = i
    for k in range(i+1, n):
        if abs(A1[k][i]) > abs(A1[max_row][i]):
            max_row = k
    A1[i], A1[max_row] = A1[max_row], A1[i]
    b1[i], b1[max_row] = b1[max_row], b1[i]

    for k in range(i+1, n):
        factor = A1[k][i] / A1[i][i]
        for j in range(i, n):
            A1[k][j] -= factor * A1[i][j]
        b1[k] -= factor * b1[i]

x1 = [0 for _ in range(n)]
for i in range(n-1, -1, -1):
    s = b1[i]
    for j in range(i+1, n):
        s -= A1[i][j] * x1[j]
    x1[i] = s / A1[i][i]

print("===== Problem 1: Gaussian Elimination Solution =====")
for i in range(n):
    print("x%d = %.6f" % (i+1, x1[i]))

# 題目二：找反矩陣
A2 = [
    [4, 1, -1, 0],
    [1, 3, -1, 0],
    [-1, -1, 6, 2],
    [0, 0, 2, 5]
]
n = 4
I = [[float(i == j) for j in range(n)] for i in range(n)]
for i in range(n):
    # Pivot
    max_row = i
    for k in range(i+1, n):
        if abs(A2[k][i]) > abs(A2[max_row][i]):
            max_row = k
    A2[i], A2[max_row] = A2[max_row], A2[i]
    I[i], I[max_row] = I[max_row], I[i]

    pivot = A2[i][i]
    for j in range(n):
        A2[i][j] /= pivot
        I[i][j] /= pivot

    for k in range(n):
        if k != i:
            factor = A2[k][i]
            for j in range(n):
                A2[k][j] -= factor * A2[i][j]
                I[k][j] -= factor * I[i][j]

print("\n===== Problem 2: Inverse of Matrix A =====")
for row in I:
    print(["%.6f" % val for val in row])

# 題目三：Crout Factorization 解三對角矩陣
A3 = [
    [3, -1, 0, 0],
    [-1, 3, -1, 0],
    [0, -1, 3, -1],
    [0, 0, -1, 3]
]
b3 = [2, 3, 4, 1]
n = 4
L = [[0]*n for _ in range(n)]
U = [[0]*n for _ in range(n)]
y = [0]*n
x3 = [0]*n

for i in range(n):
    for j in range(i+1):
        s = sum(L[i][k]*U[k][j] for k in range(j))
        L[i][j] = A3[i][j] - s if i >= j else 0
    for j in range(i, n):
        if i == j:
            U[i][j] = 1
        else:
            s = sum(L[i][k]*U[k][j] for k in range(i))
            U[i][j] = (A3[i][j] - s) / L[i][i]

for i in range(n):
    s = sum(L[i][j]*y[j] for j in range(i))
    y[i] = (b3[i] - s) / L[i][i]

for i in range(n-1, -1, -1):
    s = sum(U[i][j]*x3[j] for j in range(i+1, n))
    x3[i] = y[i] - s

print("\n===== Problem 3: Crout Factorization Solution =====")
for i in range(n):
    print("x%d = %.6f" % (i+1, x3[i]))
