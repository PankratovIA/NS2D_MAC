from math import exp, sin, cos

def taylor_u(t, X, Re):
    x, y = X
    return -exp(-2.*t/Re) * cos(x) * sin(y)
    
if __name__ == "__main__":
    print("taylor >>>")
    print("taylor <<<")
    