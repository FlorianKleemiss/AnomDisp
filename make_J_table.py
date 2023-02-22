from matrix_coefficients_v2 import *
epsilon = 1e-7

def close_to_integer(num, integer):
    is_close = False
    res = num*integer - int(num*integer)
    if (num * integer) % 1 < epsilon:
        is_close = True
    return is_close

Ws = [[],[],[],[],[],[]]
for i in range(4):
    Ws[0].append(make_matrix_W(0, i, 20))
    Ws[1].append(make_matrix_W(1, i, 20))
    Ws[2].append(make_matrix_W(2, i, 20))
    Ws[3].append(make_matrix_W(3, i, 20))

for a in range(4):
    for c in range(4):
        print("Table for J(%d,p,%d,l)"%(a,c)) 
        print("  l =        0           1           2           3           4           5           6")
        print("----------------------------------------------------------------------------------------------")
        for p in range(7):
            up = []
            down = []
            found = []
            for l in range(7):
                res = J(a,p,c,l,Ws[a][c])
                f = False
                for i in range(1,30000):
                    if close_to_integer(res,i):
                        up.append(res * i)
                        down.append(i)
                        f = True
                        break
                if f == False:
                    up.append(res)
                    down.append(-1)
                found.append(f)
            s = "p = %d|"%p
            for i in range(7):
                if found[i] == True:
                    s += "   %4d/%4d"%(up[i],down[i])
                else:
                    s += "%12.5e"%up[i]
            print(s)
        print("")
        print("")