import math
import sympy as sp
import mpmath
from sympy import re, im
from sympy.matrices import Matrix, eye

# "a" = 2*o/b
# "b" = b/2
# Force_T = F1*w**2*Y+F2*sp.I*w*PHI+F3*w**2*PHI+F4*sp.I*w*Y+F5*PHI
# Moment_T = M1*w**2*Y+M2*sp.I*w*PHI+M3*w**2*PHI+M4*sp.I*w*Y+M5*PHI
# Force_theodorsen = math.pi*ro*(b/2)**2*(-w_j**2*y+sp.I*w_j*v*phi+w_j**2*(b/2)*(1/2+o/b-1)*phi)
# Force_theodorsen+= 2*math.pi*ro*v*(b/2)*C_theod(v,w_j,b)*(sp.I*w_j*y+v*phi+sp.I*w_j*(b/2)*(1-o/b)*phi)
# Moment_theodorsen = math.pi*ro*(b/2)**2*(-w_j**2*(b/2)*(1/2+o/b-1)*y-sp.I*w_j*v*(b/2)*(1-o/b)*phi+(b/2)**2*w_j**2*(1/8+(1/2+o/b-1)**2)*phi)
# Moment_theodorsen+= 2*math.pi*ro*v*(b/2)**2*(o/b)*C_theod(v,w_j,b)*(sp.I*w_j*y+v*phi+sp.I*w_j*b*(1-o/b)*phi)

def C_theod(v,w_j,b):
    k = w_j*b/2/v
    i_k=sp.expand(sp.I*k)
    return mpmath.besselk(1, i_k)/(mpmath.besselk(0,i_k)+mpmath.besselk(1,i_k))

def theor_calc(v=150*1000, flutter_speed=0, file="result for 300_1 iter.tra"):
############ Затухание зависит от скорости узла
    v = v/3.6*1000
    Sound_speed = 340.294*1000
    file_res_300_iter = open(file,'a+')
    c1=0.09; k2=30; m=2*(10**-8); h=10; b=100; o=-30; dcy=6.283; ro=1.225e-12; l=4000; #a=(50-100*0.228);
    a=25
    A = 2*o/b
    k2 = b**2/4*c1-c1*o**2
    Iz = m*b*l*(h**2+b**2)/12
    M = m*b*l
    print(f"m = {M}; Iz = {Iz}")
    a11 = c1*l
    a12 = -o*c1*l
    a21 = -o*c1*l
    a22 = o**2*c1*l+k2*l
    ka12 = -dcy*ro*v**2/2*b*l/M
    ka22 = -a*dcy*ro*v**2/2*b*l/Iz
    a13 = dcy*ro*v/2*b*l
    a23 = dcy*ro*v/2*b*l*a
    a14 = -dcy*ro*v/2*b*l*o
    a24 = -dcy*ro*v/2*b*l*o*a
    a11/= M
    a12/= M
    a13/= M
    a14/= M
    a21/= Iz
    a22/= Iz
    a23/= Iz
    a24/= Iz
    Kk_matrix = Matrix([[a11,a12],[a21,a22]])
    Ka_matrix = Matrix([[0,ka12],[0,ka22]])
    D_matrix = Matrix([[a13,a14],[a23,a24]])
    w = sp.Symbol("w")
    matrix_eq = -w**2*eye(2)+sp.I*w*D_matrix+Kk_matrix+Ka_matrix
    eq_det = sp.Eq(matrix_eq.det(),0)
    res = sp.solve(eq_det, w)
    print(f"\n{res=}")
    if flutter_speed[0] == '>1000' and (im(res[0])<0 or im(res[1])<0 or im(res[2])<0 or im(res[3])<0):
        flutter_speed[0] = v*3.6/1000

############ Затухание зависит от скорости центра масс
    file_res_300_iter = open(file,'a+')
    a11 = c1*l
    a12 = -o*c1*l
    a21 = -o*c1*l
    a22 = o**2*c1*l+k2*l
    ka12 = -dcy*ro*v**2/2*b*l/M
    ka22 = -a*dcy*ro*v**2/2*b*l/Iz
    a13 = dcy*ro*v/2*b*l
    a23 = dcy*ro*v/2*b*l*a
    a14 = -dcy*ro*v/2*b*l*o
    a24 = -dcy*ro*v/2*b*l*o*a
    a11/= M
    a12/= M
    a13/= M
    a14/= M
    a21/= Iz
    a22/= Iz
    a23/= Iz
    a24/= Iz
    Kk_matrix = Matrix([[a11,a12],[a21,a22]])
    Ka_matrix = Matrix([[0,ka12],[0,ka22]])
    D_matrix = Matrix([[a13,0],[a23,0]])
    w = sp.Symbol("w")
    matrix_eq = -w**2*eye(2)+sp.I*w*D_matrix+Kk_matrix+Ka_matrix
    eq_det = sp.Eq(matrix_eq.det(),0)
    res = sp.solve(eq_det, w)
    print(f"\n{res=}")
    if flutter_speed[1] == '>1000' and (im(res[0])<0 or im(res[1])<0 or im(res[2])<0 or im(res[3])<0):
        flutter_speed[1] = v*3.6/1000
    

############      Theodorsen     #################
##
##    res_det = res
##    stop_flag = 1
##    ind = 0
##    while stop_flag:
##        if ind == 0:
##            w_j = res_det[0]
##            ind+= 1
##        else:
##            dist = abs(w_j_prev-res_det[0])
##            w_j = res_det[0]
##            for w_k in res_det:
##                if abs(w_j_prev-w_k)<dist:
##                    dist =abs(w_j_prev-w_k)
##                    w_j = w_k
##        w_j_prev = w_j
##        res_det_prev = res_det
####        w_j = 0.000000000001
##        print(f"{w_j=}")
##        k_j = w_j*b/2/v
####        L_z = 2*math.pi*(-(w*b/2/v)**2/2-im(C_theod(v,w_j,b))*(w*b/2/v))
####        L_z_d = 2*math.pi*sp.re(C_theod(v,w_j,b))
####        L_t = 2*math.pi*((w*b/2/v)**2*(2*o/b)/2+re(C_theod(v,w_j,b))-im(C_theod(v,w_j,b))*(w*b/2/v)*(1/2-2*o/b))
####        L_t_d = 2*math.pi*(1/2+re(C_theod(v,w_j,b))*(1/2-2*o/b)+im(C_theod(v,w_j,b))/(w*b/2/v))
####        M_z = 2*math.pi*(-(w*b/2/v)**2*(2*o/b)/2-im(C_theod(v,w_j,b))*(w*b/2/v)*(1/2+2*o/b))
####        M_z_d = 2*math.pi*(1/2+2*o/b)*re(C_theod(v,w_j,b))
####        M_t = 2*math.pi*((w*b/2/v)**2/2*(1/8+(2*o/b)**2)+re(C_theod(v,w_j,b))*(1/2+2*o/b)-im(C_theod(v,w_j,b))*(w*b/2/v)*(1/2+2*o/b)*(1/2-2*o/b))
####        M_t_d = 2*math.pi*(-1/2*(1/2-2*o/b)+re(C_theod(v,w_j,b))*(1/2+2*o/b)*(1/2-2*o/b)+im(C_theod(v,w_j,b))/(w*b/2/v)*(1/2+2*o/b))
##        L_z = 2*math.pi*(-(w_j*b/2/v)**2/2-im(C_theod(v,w_j,b))*(w_j*b/2/v))
##        L_z_d = 2*math.pi*sp.re(C_theod(v,w_j,b))
##        L_t = 2*math.pi*((w_j*b/2/v)**2*(2*o/b)/2+re(C_theod(v,w_j,b))-im(C_theod(v,w_j,b))*(w_j*b/2/v)*(1/2-2*o/b))
##        L_t_d = 2*math.pi*(1/2+re(C_theod(v,w_j,b))*(1/2-2*o/b)+im(C_theod(v,w_j,b))/(w_j*b/2/v))
##        M_z = 2*math.pi*(-(w_j*b/2/v)**2*(2*o/b)/2-im(C_theod(v,w_j,b))*(w_j*b/2/v)*(1/2+2*o/b))
##        M_z_d = 2*math.pi*(1/2+2*o/b)*re(C_theod(v,w_j,b))
##        M_t = 2*math.pi*((w_j*b/2/v)**2/2*(1/8+(2*o/b)**2)+re(C_theod(v,w_j,b))*(1/2+2*o/b)-im(C_theod(v,w_j,b))*(w_j*b/2/v)*(1/2+2*o/b)*(1/2-2*o/b))
##        M_t_d = 2*math.pi*(-1/2*(1/2-2*o/b)+re(C_theod(v,w_j,b))*(1/2+2*o/b)*(1/2-2*o/b)+im(C_theod(v,w_j,b))/(w_j*b/2/v)*(1/2+2*o/b))
##        m11 = m*b*l
##        m12 = m*b*l*o
##        m21 = m*b*l*o
##        m22 = m*b*l*o**2+Iz
##        M_matrix = Matrix([[m11,m12],[m21,m22]])
##        kk11 = c1*l
##        kk22 = k2*l
##        Kk_matrix = Matrix([[kk11,0],[0,kk22]])
####          Сейчас везде стоит w_j
####          Как превратить в квазистационарный: В строках f12 и f22 w заменить на w_j; вернуть строку w_j = 0.000000000001
##        f11 = sp.expand(ro*v**2*(b/2)*(L_z+sp.I*(w*b/2/v)*L_z_d)/(b/2))
##        f12 = sp.expand(ro*v**2*(b/2)*(L_t+sp.I*(w*b/2/v)*L_t_d))
##        f21 = sp.expand(ro*v**2*(b/2)**2*(M_z+sp.I*(w*b/2/v)*M_z_d)/(b/2))
##        f22 = sp.expand(ro*v**2*(b/2)**2*(M_t+sp.I*(w*b/2/v)*M_t_d))
##        f11*= -l
##        f12*= l
##        f21*= -l
##        f22*= l
####        print(f"{L_z = }; {L_z_d = }; {L_t = }; {L_t_d = }; {M_z = }; {M_z_d = }; {M_t = }; {M_t_d = }; ")
##        Force_matrix = Matrix([[f11,f12],[f21,f22]])
####        print(f"{Force_matrix=}")
##        Final_matrix = -w**2*M_matrix+Kk_matrix-Force_matrix
##        eq_det = sp.Eq(Final_matrix.det(),0)
##        res_det = sp.solve(eq_det, w)
##        print(f"\n{res_det=}")
##        
##        for elem in res_det:
##            close_flag = 0
##            for elem_prev in res_det_prev:
##                if abs(elem-elem_prev) <= 0.0001:
##                    close_flag = 1
##            if close_flag == 0:
##                break
##        else:
##            stop_flag = 0
##    if flutter_speed[2] == '>1000' and (im(res_det[0])<0 or im(res_det[1])<0 or im(res_det[2])<0 or im(res_det[3])<0):
##        flutter_speed[2] = v*3.6/1000
##
#################      Possio           ##################
##
##    print("###############   Possio   ##############")
##    res_det=res
##    stop_flag = 1
##    ind = 0
##    while stop_flag:
##        if ind == 0:
##            w_j = res_det[0]
##            ind+= 1
##        else:
##            dist = abs(w_j_prev-res_det[0])
##            w_j = res_det[0]
##            for w_k in res_det:
##                if abs(w_j_prev-w_k)<dist:
##                    dist =abs(w_j_prev-w_k)
##                    w_j = w_k
##        w_j_prev = w_j
##        k_j = sp.I*w_j*b/2/v
##        k = sp.I*w*b/2/v
##        res_det_prev = res_det
##        Mach = v/Sound_speed
##        Mach = 0.00000000001
##        L_h = k**2+2*k*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*(k_j**4/2+2*k_j**3*C_theod(v,w_j,b)+2*k**2*C_theod(v,w_j,b)**2)
##        L_a = k-A*k**2+C_theod(v,w_j,b)*(2+(1-2*A)*k)+Mach**2*math.log(Mach)*(k_j**3/2-k_j**4/2*A+C_theod(v,w_j,b)*(2*k**2+(1/2-2*A)*k_j**3)+C_theod(v,w_j,b)**2*(2*k+(1-2*A)*k**2))
##        M_h = A*k**2+(1+2*A)*k*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*((1+2*A)*k**2*C_theod(v,w_j,b)**2+(1/2+2*A)*k_j**3*C_theod(v,w_j,b)+k_j**4/2*A)
##        M_a = (A-1/2)*k-(1/8+A**2)*k**2+C_theod(v,w_j,b)*(1+2*A+(1/2-2*A**2)*k)+Mach**2*math.log(Mach)*(A/2*k_j**3-A**2/2*k_j**4+C_theod(v,w_j,b)**2*((1+2*A)*k+(1/2-2*A**2)*k**2)+C_theod(v,w_j,b)*((1/2+2*A)*k**2-2*A**2*k_j**3))
##
####        L_h = (sp.I*w_j*b/2/v)**2+2*(sp.I*w*b/2/v)*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*(k_j**4/2+2*k_j**3*C_theod(v,w_j,b)+2*(sp.I*w_j*b/2/v)**2*C_theod(v,w_j,b)**2)
####        L_a = (sp.I*w*b/2/v)-A*(sp.I*w_j*b/2/v)**2+C_theod(v,w_j,b)*(2+(1-2*A)*(sp.I*w*b/2/v))+Mach**2*math.log(Mach)*(k_j**3/2-k_j**4/2*A+C_theod(v,w_j,b)*(2*(sp.I*w_j*b/2/v)**2+(1/2-2*A)*k_j**3)+C_theod(v,w_j,b)**2*(2*(sp.I*w*b/2/v)+(1-2*A)*(sp.I*w_j*b/2/v)**2))
####        M_h = A*(sp.I*w_j*b/2/v)**2+(1+2*A)*(sp.I*w*b/2/v)*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*((1+2*A)*(sp.I*w_j*b/2/v)**2*C_theod(v,w_j,b)**2+(1/2+2*A)*k_j**3*C_theod(v,w_j,b)+k_j**4/2*A)
####        M_a = (A-1/2)*(sp.I*w*b/2/v)-(1/8+A**2)*(sp.I*w_j*b/2/v)**2+C_theod(v,w_j,b)*(1+2*A+(1/2-2*A**2)*(sp.I*w*b/2/v))+Mach**2*math.log(Mach)*(A/2*k_j**3-A**2/2*k_j**4+C_theod(v,w_j,b)**2*((1+2*A)*(sp.I*w*b/2/v)+(1/2-2*A**2)*(sp.I*w_j*b/2/v)**2)+C_theod(v,w_j,b)*((1/2+2*A)*(sp.I*w_j*b/2/v)**2-2*A**2*k_j**3))
##        m11 = m*b*l
##        m12 = m*b*l*o
##        m21 = m*b*l*o
##        m22 = m*b*l*o**2+Iz
##        M_matrix = Matrix([[m11,m12],[m21,m22]])
##        f11 = sp.expand(math.pi*ro*v**2*L_h)
##        f12 = sp.expand(math.pi*ro*v**2*(b/2)*L_a)
##        f21 = sp.expand(math.pi*ro*v**2*(b/2)*M_h)
##        f22 = sp.expand(math.pi*ro*v**2*(b/2)**2*M_a)
##        f11*= -l
##        f12*= l
##        f21*= -l
##        f22*= l
##        Force_matrix = Matrix([[f11,f12],[f21,f22]])
####        print(f"{Force_matrix=}")
##        Final_matrix = -w**2*M_matrix+Kk_matrix-Force_matrix
##        eq_det = sp.Eq(Final_matrix.det(),0)
##        res_det = sp.solve(eq_det, w)
##        print(f"\n{res_det=}")
##        
##        for elem in res_det:
##            close_flag = 0
##            for elem_prev in res_det_prev:
##                if abs(elem-elem_prev) <= 0.0001:
##                    close_flag = 1
##            if close_flag == 0:
##                break
##        else:
##            stop_flag = 0
##    if flutter_speed[3] == '>1000' and (im(res_det[0])<0 or im(res_det[1])<0 or im(res_det[2])<0 or im(res_det[3])<0):
##        flutter_speed[3] = v*3.6/1000
    

###############      Possio 2          ##################

    print("###############   Possio 2   ##############")
    res_det=res
    stop_flag = 1
    ind = 0
    while stop_flag:
        if ind == 0:
            w_j = res_det[0]
            ind+= 1
        else:
            dist = abs(w_j_prev-res_det[0])
            w_j = res_det[0]
            for w_k in res_det:
                if abs(w_j_prev-w_k)<dist:
                    dist =abs(w_j_prev-w_k)
                    w_j = w_k
        w_j_prev = w_j
        k_j = sp.I*w_j*b/2/v
        k = sp.I*w*b/2/v
        res_det_prev = res_det
        Mach = v/Sound_speed

        L_h = 1+2/k_j*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*(k**2/2+2*k*C_theod(v,w_j,b)+2*C_theod(v,w_j,b)**2)
        L_a = 1/2+1/k_j+C_theod(v,w_j,b)*(2/k_j**2+2/k_j)++Mach**2*math.log(Mach)*(k/2+k**2/4+C_theod(v,w_j,b)*(2+3/2*k)+C_theod(v,w_j,b)**2*(2+2/k_j))
        M_h = 1/2+Mach**2*math.log(Mach)*(k**2/4+k/2*C_theod(v,w_j,b))
        M_a = 3/8+1/k_j+Mach**2*math.log(Mach)*(k/4+k**2/8+C_theod(v,w_j,b)*(1/2+k/2))

        m11 = m*b*l
        m12 = m*b*l*o
        m21 = m*b*l*o
        m22 = m*b*l*o**2+Iz
        M_matrix = Matrix([[m11,m12],[m21,m22]])
        f11 = sp.expand(-math.pi*ro*v**2*k_j**2*L_h)
        f12 = sp.expand(-math.pi*ro*v**2*(b/2)*k_j**2*(L_a-(1/2+A)*L_h))
        f21 = sp.expand(-math.pi*ro*v**2*(b/2)*k_j**2*((M_h-(1/2+A)*L_h)))
        f22 = sp.expand(-math.pi*ro*v**2*(b/2)**2*k_j**2*(M_a*k_j-(1/2+a)*(L_a+M_h)+(1/2+a)**2*L_h))
        f11*= -l
        f12*= l
        f21*= -l
        f22*= l
        Force_matrix = Matrix([[f11,f12],[f21,f22]])
##        print(f"{Force_matrix=}")
        Final_matrix = -w**2*M_matrix+Kk_matrix-Force_matrix
        eq_det = sp.Eq(Final_matrix.det(),0)
        res_det = sp.solve(eq_det, w)
        print(f"\n{eq_det=}")
        print(f"\n{res_det=}")
        
        for elem in res_det:
            close_flag = 0
            for elem_prev in res_det_prev:
                if abs(elem-elem_prev) <= 0.0001:
                    close_flag = 1
            if close_flag == 0:
                break
        else:
            stop_flag = 0
    if flutter_speed[4] == '>1000' and (im(res_det[0])<0 or im(res_det[1])<0 or im(res_det[2])<0 or im(res_det[3])<0):
        flutter_speed[4] = v*3.6/1000
    return flutter_speed
if __name__ == "__main__":
    flutter_speed = ['>1000', '>1000', '>1000', '>1000', '>1000']
    for i in range(30, 500, 20):
##        file_res_300_iter = open("result for 300_1 iter.tra",'a+')
####        file_res_300_iter.write(str(i))
##        file_res_300_iter.close()
        flutter_speed = theor_calc(i, flutter_speed)
        print(f"\n\nspeed = {i}, flutter = {flutter_speed}\n")
