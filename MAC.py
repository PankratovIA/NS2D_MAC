from math import exp, sin, cos, pi
import numpy as np

Re = 100.
h = 0.2
tau = 0.05 * h

a = b = 1
cnt = int(a / h)
size = cnt * cnt

def taylor_u(t, X):
    x, y = X
    return -exp(-2.*t/Re) * cos(x) * sin(y)
    
def taylor_v(t, X):
    x, y = X
    return exp(-2.*t/Re) * sin(x) * cos(y)

def taylor_p(t, X):
    x, y = X
    # print(x, y)
    return -exp(-4.*t/Re) * (cos(2.*x) + cos(2.*y))/4
    
def num(row, col, cnt):
    return row * cnt + col

def F(j, k, u, v):
    # F_{j+1/2, k}
    # u[j][k] = u_{j+1/2, k}; v[j][k] = v_{j, k+1/2}; p[j][k] = p_{j, k}
    ans = (u[j+1][k] - 2 * u[j][k] + u[j-1][k]) / (Re * (h ** 2)) # dx**2
    ans += (u[j][k-1] - 2 * u[j][k] + u[j][k+1]) / (Re * (h ** 2)) # dy**2
    
    # u^2
    u_sqr = ((u[j][k] + u[j+1][k]) ** 2) / 4.0
    u_sqr -= ((u[j-1][k] + u[j][k]) ** 2) / 4.0
    ans -= u_sqr / h
    
    # uv
    uv = 0.25 * (u[j][k] + u[j][k+1]) * (v[j+1][k] + v[j][k])
    uv -= 0.25 * (u[j][k-1] + u[j][k]) * (v[j+1][k-1] + v[j][k-1])
    ans -= uv / h
    
    ans *= tau
    ans += u[j][k]
    return ans
    
def G(j, k, u, v):
    # G_{j, k+1/2}
    # u[j][k] = u_{j+1/2, k}; v[j][k] = v_{j, k+1/2}; p[j][k] = p_{j, k}
    ans = (v[j+1][k] - 2 * v[j][k] + v[j-1][k]) / (Re * (h ** 2)) # dx**2
    ans += (v[j][k+1] - 2 * v[j][k] + v[j][k-1]) / (Re * (h ** 2)) # dy**2

    # uv
    uv = 0.25 * (u[j][k] + u[j][k+1]) * (v[j+1][k] + v[j][k])
    uv -= 0.25 * (u[j-1][k] + u[j-1][k+1]) * (v[j][k] + v[j-1][k])
    ans -= uv / h

    # v^2
    v_sqr = ((v[j][k] + v[j][k+1]) ** 2) / 4.0
    v_sqr -= ((v[j][k-1] + v[j][k]) ** 2) / 4.0
    ans -= v_sqr / h
  
    ans *= tau
    ans += v[j][k]
    return ans
    
def genBC(t, x, y, func):
    A = np.zeros((size, size))
    b = np.zeros((size, 1))
    for cur in range(cnt):
        col, row = cur, cur
        for tmp in range(2):
            #bottom
            idx = num(0 + tmp, col, cnt)
            A[idx][idx] = 1.0
            b[idx] = func(t, (x[col], y[0]))
            #top
            idx = num(cnt - 1 - tmp, col, cnt)
            A[idx][idx] = 1.0
            b[idx] = func(t, (x[col], y[cnt - 1]))
            #left
            idx = num(row, 0 + tmp, cnt)
            A[idx][idx] = 1.0
            b[idx] = func(t, (x[0], y[row]))
            #right
            idx = num(row, cnt - 1 - tmp, cnt)
            A[idx][idx] = 1.0
            b[idx] = func(t, (x[cnt - 1], y[row]))
    return A, b
    
def genMatrP(t, x, y, u, v):
    A, b = genBC(t, x, y, taylor_p)

    for j in range(2, cnt-2):
        for k in range(2, cnt-2):
            row = num(j, k, cnt)   
            print("eq p", j, k, row)
            # dx
            A[row][row] += -2.0
            col = num(j-1, k, cnt)
            A[row][col] += 1.0
            col = num(j+1, k, cnt)
            A[row][col] += 1.0
            # dy
            A[row][row] += -2.0
            col = num(j, k-1, cnt)
            A[row][col] += 1.0
            col = num(j, k+1, cnt)
            A[row][col] += 1.0
            
            b[row] = F(j, k, u, v) - F(j-1, k, u, v)
            b[row] += G(j, k, u, v) - G(j, k-1, u, v)
            b[row] *= h/tau
    # print(A)
    # print(b)
    
    print(A[12])
    print(b[12])
    return A, b
    
if __name__ == "__main__":
    print("taylor >>>")
    #print(taylor_v(0, (h, h), Re))
    
    p_x = p_y = [(j + 0.5) * h for j in range(cnt)]
    print("p_x =", p_x)
    print("p_y =", p_y)
    
    print("p, u, v t=0 >>>")
    t = 0.
    p = [[taylor_p(t, (p_x[j], p_y[k]))for j in range(cnt)] for k in range(cnt)]
    
    for row in p:
        print(" ".join(["{0:7.4f}".format(cur)  for cur in row]))
    print(taylor_p(0, (h/2, h/2)))
    print(taylor_p(0, (3.*h/2, 5.*h/2)), '\n')
    
    u_x = [ (j + 1) * h for j in range(cnt)]
    u_y = [(j + 0.5) * h for j in range(cnt)] # as p_y
    print("u_x =", u_x)
    print("u_y =", u_y)
    
    t = 0.
    u = [[taylor_u(t, (u_x[j], u_y[k]))for j in range(cnt)] for k in range(cnt)]
    for row in u:
        print(" ".join(["{0:7.4f}".format(cur)  for cur in row]))
    print(taylor_u(0, (h, h/2)))
    print(taylor_u(0, (3.*h, 5.*h +h/2)), '\n')
    
    v_x = [ (j + 0.5) * h for j in range(cnt)] # as p_x
    v_y = [(j + 1) * h for j in range(cnt)] 
    print("v_x =", v_x)
    print("v_y =", v_y)
    
    t = 0.
    v = [[taylor_v(t, (v_x[j], v_y[k]))for j in range(cnt)] for k in range(cnt)]
    for row in v:
        print(" ".join(["{0:7.4f}".format(cur)  for cur in row]))
    print(taylor_v(0, (h/2, h)))
    print(taylor_v(0, (3.*h + h/2, 5.*h)), '\n')
    print("p, u, v t=0 <<<")
    
    Ap, bp = genMatrP(t, p_x, p_y, u, v)
    print("p det = ", np.linalg.det(Ap))
    p1 = np.linalg.solve(Ap, bp)
    print(p1, p1[12])
    
    
    print("taylor <<<")

