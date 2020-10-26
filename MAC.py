from math import exp, sin, cos, pi
import numpy as np

Re = 100.
h = 0.1
tau = 0.05 * h

a = b = 1
cnt = int(a / h)

def taylor_u(t, X, Re):
    x, y = X
    return -exp(-2.*t/Re) * cos(x) * sin(y)
    
def taylor_v(t, X, Re):
    x, y = X
    return exp(-2.*t/Re) * sin(x) * cos(y)

def taylor_p(t, X, Re):
    x, y = X
    # print(x, y)
    return -exp(-4.*t/Re) * (cos(2.*x) + cos(2.*y))/4
    
if __name__ == "__main__":
    print("taylor >>>")
    #print(taylor_v(0, (h, h), Re))
    
    p_x = p_y = [(j + 0.5) * h for j in range(cnt)]
    print("p_x =", p_x)
    print("p_y =", p_y)
    
    t = 0.
    p = [[taylor_p(t, (p_x[j], p_y[k]), Re)for j in range(cnt)] for k in range(cnt)]
    
    for row in p:
        print(" ".join(["{0:7.4f}".format(cur)  for cur in row]))
    print(taylor_p(0, (h/2, h/2), Re))
    print(taylor_p(0, (3.*h/2, 5.*h/2), Re), '\n')
    
    u_x = [ (j + 1) * h for j in range(cnt)]
    u_y = [(j + 0.5) * h for j in range(cnt)] # as p_y
    print("u_x =", u_x)
    print("u_y =", u_y)
    
    t = 0.
    u = [[taylor_u(t, (u_x[j], u_y[k]), Re)for j in range(cnt)] for k in range(cnt)]
    for row in u:
        print(" ".join(["{0:7.4f}".format(cur)  for cur in row]))
    print(taylor_u(0, (h, h/2), Re))
    print(taylor_u(0, (3.*h, 5.*h +h/2), Re), '\n')
    
    v_x = [ (j + 0.5) * h for j in range(cnt)] # as p_x
    v_y = [(j + 1) * h for j in range(cnt)] 
    print("v_x =", v_x)
    print("v_y =", v_y)
    
    t = 0.
    v = [[taylor_v(t, (v_x[j], v_y[k]), Re)for j in range(cnt)] for k in range(cnt)]
    for row in v:
        print(" ".join(["{0:7.4f}".format(cur)  for cur in row]))
    print(taylor_v(0, (h/2, h), Re))
    print(taylor_v(0, (3.*h + h/2, 5.*h), Re), '\n')
    
    print("taylor <<<")