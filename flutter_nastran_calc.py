import math
import sympy as sp
import mpmath
import subprocess
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

def Solve_nastran(Nastran_matrix_1_deg,Nastran_matrix_0_deg):
    file = open("DMIG_file.pch", 'w')
    write_DMIG_KGG()
    write_DMIG(Nastran_matrix_0_deg[0,0], 0, 2, 2)
    write_DMIG(Nastran_matrix_0_deg[0,1], 0, 2, 6)
    write_DMIG(Nastran_matrix_0_deg[1,0], 0, 6, 2)
    write_DMIG(Nastran_matrix_0_deg[1,1], 0, 6, 6)
    write_DMIG_DGG()
    write_DMIG(Nastran_matrix_1_deg[0,0], 1, 2, 2)
    write_DMIG(Nastran_matrix_1_deg[0,1], 1, 2, 6)
    write_DMIG(Nastran_matrix_1_deg[1,0], 1, 6, 2)
    write_DMIG(Nastran_matrix_1_deg[1,1], 1, 6, 6)
    nastran_commands = [r"C:\MSC.Software\MSC_Nastran\20160\bin\nastran.exe", "wing_conm_sol_107.bdf", r"currdir=" + r"C:\Users\New\Desktop\programm\test2" ]
    result = subprocess.run(nastran_commands, shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    res_sol_107 = open("wing_conm_sol_107.f06",'r')
    res = []
    start_of_res = 0
    for line in res_sol_107:
        if start_of_res == 0:
            if "C O M P L E X   E I G E N V A L U E   S U M M A R Y" in line:
                start_of_res = 1
##                print(line)
        elif start_of_res <=2:
            start_of_res += 1
        elif start_of_res <=6:
            start_of_res +=1
##            print(line)
            res.append(float(line.split()[2])+sp.I*float(line.split()[3]))
##            print(res)
##    print(f"{res=}")
    res_sol_107.close()
    return res
def write_DMIG(value, degree, row, column):
    file = open("DMIG_file.pch", 'a')
    file.write(f"DMIG*   {'KGG400S' if degree == 0 else 'DGG400S'}         1               {column}\n")
    file.write(f"*       1               {row}               {float(re(value)):<16.7f}{float(im(value)):<16.7f}\n")
def write_DMIG_KGG():
    file = open("DMIG_file.pch", 'a')
    file.write(f"DMIG    KGG400S 0       1       3       0\n")
def write_DMIG_DGG():
    file = open("DMIG_file.pch", 'a')
    file.write(f"DMIG    DGG400S 0       1       3       0\n")
def create_bdf_file(c1):
    k2=b**2/4*c1-c1*o**2
    bdf = open("wing_conm_sol_107.bdf", 'w')
    bdf.write(f"""$$$$$$$$$$$$$$$$$$$$$$$$$$$$ The NASTRAN Statement $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
NASTRAN SYSTEM(316)=19
NASTRAN SYSTEM(444)=9
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$ File Management Statements $$$$$$$$$$$$$$$$$$$$$$$$$$
$ASSIGN INPUTT4='matrixd_case_50005.1.op4',OLD,UNIT=23,FORMATTED
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$ Executive Control Statements $$$$$$$$$$$$$$$$$$$$$$$$$
SOL 107
$ DMAP
DIAG = 4,8,14,15
COMPILE SEDCEIG
alter 22
$INPUTT4   /KGGSQ,MGGSQ,BDAM,,/3/23/-1// $
$UPARTN   USET,KGGSQ/KAASQ,,,/'G'/'A'// $
$UPARTN   USET,MGGSQ/MAASQ,,,/'G'/'A'// $
$UMERGE1  USET,BDAM,,,/BDAG/'G'/'A'
$EQUIVX   KGGSQ/KGG/-1 $
$EQUIVX   MGGSQ/MGG/-1 $
$EQUIVX   KAASQ/KAA/-1 $
$EQUIVX   MAASQ/MAA/-1 $
$EQUIVX   BDAM/BAA/-1 $
$EQUIVX   BDAG/BGG/-1 $
$
COMPILE CEIGRS
alter 90

MATPRN    CKDXX,CBDXX,CMDXX / / $
CEND
K2PP = KGG400S
B2PP = DGG400S
SUBCASE 1
 SUBTITLE=aero
 DISPLACEMENT(SORT1,PHASE)=ALL
 SPCFORCES(SORT1,REAL)=ALL
 force(plot)=all
 CMETHOD = 10
BEGIN BULK
PARAM,grdpnt,0
PARAM,POST,0
PARAM,COUPMASS,1
PARAM,PRTMAXIM,YES
PARAM,LGDISP,1
GRID,1000,,-15.0,0.0,0.0,,123456
include 'DMIG_file.pch'
MDLPRM,OFFDEF,LROFF
EIGC, 10, CLAN, MAX, , , , , ,+ELAN1
+ELAN1, 0.0,1.0, , , , , 20,, +ELAN2



$blade start
CORD2R,1,0,0.,0.,0.,0.,0.,1580.,
 ,1.0,0.0,790.
$blade end
PBUSH*   1               K               0.             {c1*4000}
*                                                       {k2*4000}
CBUSH*   1               1               1
*                                                        0
$blade mass start
CONM2    4       2       0       0.008   0.      0.      0.
         0.      0.      0       0.      0.      6.733333
$blade mass end
RBE2*    3               1               123456          2
*
GRID*    1               1               0.              0.
*        0.              1
GRID*    2               1               -20.            0.
*        0.              1
SPC1*    1               1345            1
*
CORD2R,1,0,0.,0.,0.,0.,0.,1580.,
 ,1.0,0.0,790.
ENDDATA f0ad8b32""")
##########      Theodorsen     #################
def get_Theodorsen(v, c1, m, h, b, o, dcy, ro, l, a, Sound_speed=340.294*1000):
    
    v = v/3.6*1000
    k2 = b**2/4*c1-c1*o**2
##    print(f"{c1*4000 = }, {k2*4000 = }")
    Iz = m*b*l*(h**2+b**2)/12
    M = m*b*l
    A = 2*o/b
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
##    print(f"{res = }")
    Force_aero = dcy*ro*v**2/2*b*l
    Moment_aero = Force_aero*a
    m11 = m*b*l
    m12 = m*b*l*o
    m21 = m*b*l*o
    m22 = m*b*l*o**2+Iz
    M_matrix = Matrix([[m11,m12],[m21,m22]])
    kk11 = c1*l
    kk22 = k2*l
    ka12 = -Force_aero
    ka22 = -Force_aero*o-Moment_aero
    Kk_matrix = Matrix([[kk11,0],[0,kk22]])
    Ka_matrix = Matrix([[0,ka12],[0,ka22]])
    d11 = Force_aero/v
    d21 = Force_aero/v*o+Moment_aero/v
    D_matrix = Matrix([[d11,0],[d21,0]])
    Final_matrix = -w**2*M_matrix+w*sp.I*D_matrix+Kk_matrix+Ka_matrix
    Det=Final_matrix.det()
    equation_det = sp.Eq(Det,0)
    res_det = sp.solve(equation_det, w)
    res_det = [sp.expand(res_det[i]) for i in range(len(res_det))]
    print(f"{res_det = }")



    flutter_flag = 0
    res_det = res    
    w_j = res_det[0]
    k_j = w_j*b/2/v
    k = w*b/2/v
    L_z = 2*math.pi*(-k_j**2/2-im(C_theod(v,w_j,b))*k_j)
    L_z_d = 2*math.pi*sp.re(C_theod(v,w_j,b))
    L_t = 2*math.pi*(k_j**2*(2*o/b)/2+re(C_theod(v,w_j,b))-im(C_theod(v,w_j,b))*k_j*(1/2-2*o/b))
    L_t_d = 2*math.pi*(1/2+re(C_theod(v,w_j,b))*(1/2-2*o/b)+im(C_theod(v,w_j,b))/k_j)
    M_z = 2*math.pi*(-k_j**2*(2*o/b)/2-im(C_theod(v,w_j,b))*k_j*(1/2+2*o/b))
    M_z_d = 2*math.pi*(1/2+2*o/b)*re(C_theod(v,w_j,b))
    M_t = 2*math.pi*(k_j**2/2*(1/8+(2*o/b)**2)+re(C_theod(v,w_j,b))*(1/2+2*o/b)-im(C_theod(v,w_j,b))*k_j*(1/2+2*o/b)*(1/2-2*o/b))
    M_t_d = 2*math.pi*(-1/2*(1/2-2*o/b)+re(C_theod(v,w_j,b))*(1/2+2*o/b)*(1/2-2*o/b)+im(C_theod(v,w_j,b))/k_j*(1/2+2*o/b))
    m11 = m*b*l
    m12 = m*b*l*o
    m21 = m*b*l*o
    m22 = m*b*l*o**2+Iz
    M_matrix = Matrix([[m11,m12],[m21,m22]])
    kk11 = c1*l
    kk12 = 0
    kk22 = k2*l
    kk21 = 0
    Kk_matrix = Matrix([[kk11,0],[0,kk22]])
    f11 = -sp.expand(ro*v**2*(b/2)*(L_z+sp.I*k*L_z_d)/(b/2))
##    print(f"{f11 = }")
    f12 = -sp.expand(ro*v**2*(b/2)*(L_t+sp.I*k*L_t_d))
    f21 = -sp.expand(ro*v**2*(b/2)**2*(M_z+sp.I*k*M_z_d)/(b/2))
    f22 = -sp.expand(ro*v**2*(b/2)**2*(M_t+sp.I*k*M_t_d))
    f11*= -l
    f12*= l
    f21*= -l
    f22*= l
    Nastran_matrix_0_deg_11 = sp.expand(sp.poly(f11).coeffs()[1 if len(sp.poly(f11).coeffs()) == 2 else 0])
    Nastran_matrix_0_deg_12 = sp.expand(sp.poly(f12).coeffs()[1 if len(sp.poly(f12).coeffs()) == 2 else 0])
    Nastran_matrix_0_deg_21 = sp.expand(sp.poly(f21).coeffs()[1 if len(sp.poly(f21).coeffs()) == 2 else 0])
    Nastran_matrix_0_deg_22 = sp.expand(sp.poly(f22).coeffs()[1 if len(sp.poly(f22).coeffs()) == 2 else 0])
    Nastran_matrix_1_deg_11 = sp.expand((sp.poly(f11).coeffs()[0] if len(sp.poly(f11).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_1_deg_12 = sp.expand((sp.poly(f12).coeffs()[0] if len(sp.poly(f12).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_1_deg_21 = sp.expand((sp.poly(f21).coeffs()[0] if len(sp.poly(f21).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_1_deg_22 = sp.expand((sp.poly(f22).coeffs()[0] if len(sp.poly(f22).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_0_deg = Matrix([[Nastran_matrix_0_deg_11, Nastran_matrix_0_deg_12],[Nastran_matrix_0_deg_21, Nastran_matrix_0_deg_22]])
    Nastran_matrix_1_deg = Matrix([[Nastran_matrix_1_deg_11, Nastran_matrix_1_deg_12],[Nastran_matrix_1_deg_21, Nastran_matrix_1_deg_22]])
##    print(f"{Nastran_matrix_0_deg = }")
##    print(f"{Nastran_matrix_1_deg = }")
    Force_matrix = Matrix([[f11,f12],[f21,f22]])
##    print(f"{w*Nastran_matrix_1_deg=}\n")
    Final_matrix = -w**2*M_matrix+Kk_matrix+Force_matrix
##    print(f"{Kk_matrix+Force_matrix=}")
##    print(f"{w*sp.I*Nastran_matrix_1_deg+Nastran_matrix_0_deg+Kk_matrix=}")
    Final_matrix_nastran = w**2*M_matrix+w*Nastran_matrix_1_deg+Nastran_matrix_0_deg+Kk_matrix
##    print(f"{Final_matrix_nastran=}")
##    print(f"{Final_matrix=}")
##    eq_nastran = sp.Eq(Final_matrix_nastran.det(),0)
##    res_nastran = sp.solve(eq_nastran, w)
##    print(f"{res_nastran=}")
##    eq_det = sp.Eq(Final_matrix.det(),0)
    res_det = Solve_nastran(Nastran_matrix_1_deg,Nastran_matrix_0_deg)
    res_det = [sp.expand(i/sp.I) for i in res_det]
##    res_matrix = sp.solve(eq_det)
    print(f"{res_det=}")
##    print(f"{res_matrix=}")
    res_after_first_iteration = res_det[:]
##    print(f"{res_det = }")
##    print(f"{res_after_first_iteration=}")
##    print(f"{res_det=}")
    two_closest = []
    w_j_final = []
    closest = res_after_first_iteration[0]
    for r in res_after_first_iteration:
        if abs(r-w_j)<abs(closest-w_j):
            closest = r
    two_closest.append(closest)
    res_after_first_iteration.remove(closest)
    closest = res_after_first_iteration[0]
    for r in res_after_first_iteration:
        if abs(r-w_j)<abs(closest-w_j):
            closest = r
    two_closest.append(closest)
    res_after_first_iteration.remove(closest)
##    print(f"{res_det=}")      
    for r in two_closest:
        stop_flag = 1
        ind = 0
        while True:
##            print(f"{res_det=}")
            if ind == 0:
                ind+= 1
                w_j = r
            else:
                ind+= 1
                w_j_prev = w_j
                res_det_prev = res_det
##                print(f"{res_det=}")

                w_j = res_det[0]
                dist = abs(w_j_prev-w_j)
                for w_k in res_det:
                    if abs(w_j_prev-w_k)<dist:
                        dist = abs(w_j_prev-w_k)
                        w_j = w_k
##                    print(f"Dist from {w_j_prev} to {w_k}  is  {abs(w_j_prev-w_k)}")
    ##        w_j = 0.000000000001
##            print(f"{w_j=}")
            k_j = w_j*b/2/v
            k = w*b/2/v
            L_z = 2*math.pi*(-k_j**2/2-im(C_theod(v,w_j,b))*k)
            L_z_d = 2*math.pi*sp.re(C_theod(v,w_j,b))
            L_t = 2*math.pi*(k_j**2*(2*o/b)/2+re(C_theod(v,w_j,b))-im(C_theod(v,w_j,b))*k*(1/2-2*o/b))
            L_t_d = 2*math.pi*(1/2+re(C_theod(v,w_j,b))*(1/2-2*o/b)+im(C_theod(v,w_j,b))/k)
            M_z = 2*math.pi*(-k_j**2*(2*o/b)/2-im(C_theod(v,w_j,b))*k*(1/2+2*o/b))
            M_z_d = 2*math.pi*(1/2+2*o/b)*re(C_theod(v,w_j,b))
            M_t = 2*math.pi*(k_j**2/2*(1/8+(2*o/b)**2)+re(C_theod(v,w_j,b))*(1/2+2*o/b)-im(C_theod(v,w_j,b))*k*(1/2+2*o/b)*(1/2-2*o/b))
            M_t_d = 2*math.pi*(-1/2*(1/2-2*o/b)+re(C_theod(v,w_j,b))*(1/2+2*o/b)*(1/2-2*o/b)+im(C_theod(v,w_j,b))/k*(1/2+2*o/b))
            m11 = m*b*l
            m12 = m*b*l*o
            m21 = m*b*l*o
            m22 = m*b*l*o**2+Iz
            M_matrix = Matrix([[m11,m12],[m21,m22]])
            kk11 = c1*l
            kk22 = k2*l
            Kk_matrix = Matrix([[kk11,0],[0,kk22]])
    ##          Сейчас везде стоит w_j
    ##          Как превратить в квазистационарный: В строках f12 и f22 w заменить на w_j; вернуть строку w_j = 0.000000000001
            f11 = -sp.expand(ro*v**2*(b/2)*(L_z+sp.I*k*L_z_d)/(b/2))
            f12 = -sp.expand(ro*v**2*(b/2)*(L_t+sp.I*k*L_t_d))
            f21 = -sp.expand(ro*v**2*(b/2)**2*(M_z+sp.I*k*M_z_d)/(b/2))
            f22 = -sp.expand(ro*v**2*(b/2)**2*(M_t+sp.I*k*M_t_d))
            f11*= -l
            f12*= l
            f21*= -l
            f22*= l
            Nastran_matrix_0_deg_11 = sp.expand(sp.poly(f11).coeffs()[1 if len(sp.poly(f11).coeffs()) == 2 else 0])
            Nastran_matrix_0_deg_12 = sp.expand(sp.poly(f12).coeffs()[1 if len(sp.poly(f12).coeffs()) == 2 else 0])
            Nastran_matrix_0_deg_21 = sp.expand(sp.poly(f21).coeffs()[1 if len(sp.poly(f21).coeffs()) == 2 else 0])
            Nastran_matrix_0_deg_22 = sp.expand(sp.poly(f22).coeffs()[1 if len(sp.poly(f22).coeffs()) == 2 else 0])
            Nastran_matrix_1_deg_11 = sp.expand((sp.poly(f11).coeffs()[0] if len(sp.poly(f11).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_1_deg_12 = sp.expand((sp.poly(f12).coeffs()[0] if len(sp.poly(f12).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_1_deg_21 = sp.expand((sp.poly(f21).coeffs()[0] if len(sp.poly(f21).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_1_deg_22 = sp.expand((sp.poly(f22).coeffs()[0] if len(sp.poly(f22).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_0_deg = Matrix([[Nastran_matrix_0_deg_11, Nastran_matrix_0_deg_12],[Nastran_matrix_0_deg_21, Nastran_matrix_0_deg_22]])
            Nastran_matrix_1_deg = Matrix([[Nastran_matrix_1_deg_11, Nastran_matrix_1_deg_12],[Nastran_matrix_1_deg_21, Nastran_matrix_1_deg_22]])
    ##        print(f"{L_z = }; {L_z_d = }; {L_t = }; {L_t_d = }; {M_z = }; {M_z_d = }; {M_t = }; {M_t_d = }; ")
            Force_matrix = Matrix([[f11,f12],[f21,f22]])
    ##        print(f"{Force_matrix=}")
            Final_matrix = -w**2*M_matrix+Kk_matrix+Force_matrix
            eq_det = sp.Eq(Final_matrix.det(),0)
##            print(f"{eq_det=}")
##            res_det = sp.solve(eq_det)
##            print(f"{res_det=}")
            res_det = Solve_nastran(Nastran_matrix_1_deg,Nastran_matrix_0_deg)
            res_det = [i/sp.I for i in res_det]
##            print(f"{res_det=}")
    ##        print(f"\n{eq_det=}")
    ##        print(f"\n{res_det=}\n\n")
##            print(f"{w_j=}")
            if ind > 1 and abs(w_j-w_j_prev)<=0.001:
##                print(f" ")
                w_j_final.append(w_j)
                if im(w_j)<0:
                    flutter_flag = 1
                break
            if ind == 20:
                flutter_flag = 1
                w_j_final.append(w_j)
                print(f"Итерационный процесс не сходится. Принято, что взята скорость выше критической. Результат двух последних итераций:\n{w_j_prev}\n{w_j}  (Флаттер - отриц. мнимая часть)")
                break
    if len(w_j_final)<2:
        print("Схождение к одному корню")
    print(w_j_final)
    print(f"{f21 =};   {f22 =}")
    exit()
    return flutter_flag


#################      Possio    ################################

##    print("###############   Possio   ##############")
def get_Possio(v, c1, m, h, b, o, dcy, ro, l, a, Sound_speed=340.294*1000):
    v = v/3.6*1000
    k2 = b**2/4*c1-c1*o**2
    Iz = m*b*l*(h**2+b**2)/12
    M = m*b*l
    A = 2*o/b
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
    res = [i*sp.I for i in res]
##    print(f"{res=}")
    flutter_flag = 0
    res_det = res    
    w_j = res_det[0]
##    print(f"{w_j=}")
    k_j = sp.I*w_j*b/2/v
    k = sp.I*w*b/2/v  
    Mach = v/Sound_speed
##    Mach = 0.00000000001
    L_h = k_j**2+2*k*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*(k_j**4/2+2*k_j**3*C_theod(v,w_j,b)+2*k_j**2*C_theod(v,w_j,b)**2)
    L_a = k-A*k_j**2+C_theod(v,w_j,b)*(2+(1-2*A)*k)+Mach**2*math.log(Mach)*(k_j**3/2-k_j**4/2*A+C_theod(v,w_j,b)*(2*k_j**2+(1/2-2*A)*k_j**3)+C_theod(v,w_j,b)**2*(2*k_j+(1-2*A)*k_j**2))
    M_h = A*k_j**2+(1+2*A)*k*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*((1+2*A)*k_j**2*C_theod(v,w_j,b)**2+(1/2+2*A)*k_j**3*C_theod(v,w_j,b)+k_j**4/2*A)
    M_a = (A-1/2)*k-(1/8+A**2)*k_j**2+C_theod(v,w_j,b)*(1+2*A+(1/2-2*A**2)*k)+Mach**2*math.log(Mach)*(A/2*k_j**3-A**2/2*k_j**4+C_theod(v,w_j,b)**2*((1+2*A)*k_j+(1/2-2*A**2)*k_j**2)+C_theod(v,w_j,b)*((1/2+2*A)*k_j**2-2*A**2*k_j**3))
    m11 = m*b*l
    m12 = m*b*l*o
    m21 = m*b*l*o
    m22 = m*b*l*o**2+Iz
    M_matrix = Matrix([[m11,m12],[m21,m22]])
    kk11 = c1*l
    kk22 = k2*l
    Kk_matrix = Matrix([[kk11,0],[0,kk22]])
    f11 = -sp.expand(math.pi*ro*v**2*L_h)
    f12 = -sp.expand(math.pi*ro*v**2*(b/2)*L_a)
    f21 = -sp.expand(math.pi*ro*v**2*(b/2)*M_h)
    f22 = -sp.expand(math.pi*ro*v**2*(b/2)**2*M_a)
    f11*= -l
    f12*= l
    f21*= -l
    f22*= l
    Nastran_matrix_0_deg_11 = sp.expand(sp.poly(f11).coeffs()[1 if len(sp.poly(f11).coeffs()) == 2 else 0])
    Nastran_matrix_0_deg_12 = sp.expand(sp.poly(f12).coeffs()[1 if len(sp.poly(f12).coeffs()) == 2 else 0])
    Nastran_matrix_0_deg_21 = sp.expand(sp.poly(f21).coeffs()[1 if len(sp.poly(f21).coeffs()) == 2 else 0])
    Nastran_matrix_0_deg_22 = sp.expand(sp.poly(f22).coeffs()[1 if len(sp.poly(f22).coeffs()) == 2 else 0])
    Nastran_matrix_1_deg_11 = sp.expand((sp.poly(f11).coeffs()[0] if len(sp.poly(f11).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_1_deg_12 = sp.expand((sp.poly(f12).coeffs()[0] if len(sp.poly(f12).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_1_deg_21 = sp.expand((sp.poly(f21).coeffs()[0] if len(sp.poly(f21).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_1_deg_22 = sp.expand((sp.poly(f22).coeffs()[0] if len(sp.poly(f22).coeffs()) == 2 else 0))/sp.I
    Nastran_matrix_0_deg = Matrix([[Nastran_matrix_0_deg_11, Nastran_matrix_0_deg_12],[Nastran_matrix_0_deg_21, Nastran_matrix_0_deg_22]])
    Nastran_matrix_1_deg = Matrix([[Nastran_matrix_1_deg_11, Nastran_matrix_1_deg_12],[Nastran_matrix_1_deg_21, Nastran_matrix_1_deg_22]])
    Force_matrix = Matrix([[f11,f12],[f21,f22]])
    Final_matrix = -w**2*M_matrix+Kk_matrix+Force_matrix
##    print(f"{Kk_matrix=}")
    eq_det = sp.Eq(Final_matrix.det(),0)
    res_det = Solve_nastran(Nastran_matrix_1_deg,Nastran_matrix_0_deg)
    res_det = [i/sp.I for i in res_det]
    res_after_first_iteration = res_det
##    print(f"{res_after_first_iteration=}")
    two_closest = []
    w_j_final = []
    closest = res_after_first_iteration[0]
    for r in res_after_first_iteration:
        if abs(r-w_j)<abs(closest-w_j):
            closest = r
    two_closest.append(closest)
    res_after_first_iteration.remove(closest)
    closest = res_after_first_iteration[0]
    for r in res_after_first_iteration:
        if abs(r-w_j)<abs(closest-w_j):
            closest = r
    two_closest.append(closest)
    res_after_first_iteration.remove(closest)
          
    for r in two_closest:
        stop_flag = 1
        ind = 0
        while True:
            if ind == 0:
                ind+= 1
                w_j = r
            else:
                ind+= 1
                w_j_prev = w_j
                res_det_prev = res_det
                w_j = res_det[0]
                dist = abs(w_j_prev-w_j)
                for w_k in res_det:
                    if abs(w_j_prev-w_k)<dist:
                        dist = abs(w_j_prev-w_k)
                        w_j = w_k
##                    print(f"Dist from {w_j_prev} to {w_k}  is  {abs(w_j_prev-w_k)}")
    ##        w_j = 0.000000000001
##            print(f"{w_j=}")
            k_j = sp.I*w_j*b/2/v
            k = sp.I*w*b/2/v  
            Mach = v/Sound_speed
##            Mach = 0.00000000001
            L_h = k_j**2+2*k*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*(k_j**4/2+2*k_j**3*C_theod(v,w_j,b)+2*k_j**2*C_theod(v,w_j,b)**2)
            L_a = k-A*k_j**2+C_theod(v,w_j,b)*(2+(1-2*A)*k)+Mach**2*math.log(Mach)*(k_j**3/2-k_j**4/2*A+C_theod(v,w_j,b)*(2*k_j**2+(1/2-2*A)*k_j**3)+C_theod(v,w_j,b)**2*(2*k_j+(1-2*A)*k_j**2))
            M_h = A*k_j**2+(1+2*A)*k*C_theod(v,w_j,b)+Mach**2*math.log(Mach)*((1+2*A)*k_j**2*C_theod(v,w_j,b)**2+(1/2+2*A)*k_j**3*C_theod(v,w_j,b)+k_j**4/2*A)
            M_a = (A-1/2)*k-(1/8+A**2)*k_j**2+C_theod(v,w_j,b)*(1+2*A+(1/2-2*A**2)*k)+Mach**2*math.log(Mach)*(A/2*k_j**3-A**2/2*k_j**4+C_theod(v,w_j,b)**2*((1+2*A)*k_j+(1/2-2*A**2)*k_j**2)+C_theod(v,w_j,b)*((1/2+2*A)*k_j**2-2*A**2*k_j**3))
            m11 = m*b*l
            m12 = m*b*l*o
            m21 = m*b*l*o
            m22 = m*b*l*o**2+Iz
            M_matrix = Matrix([[m11,m12],[m21,m22]])
            kk11 = c1*l
            kk22 = k2*l
            Kk_matrix = Matrix([[kk11,0],[0,kk22]])
            f11 = -sp.expand(math.pi*ro*v**2*L_h)
            f12 = -sp.expand(math.pi*ro*v**2*(b/2)*L_a)
            f21 = -sp.expand(math.pi*ro*v**2*(b/2)*M_h)
            f22 = -sp.expand(math.pi*ro*v**2*(b/2)**2*M_a)
            f11*= -l
            f12*= l
            f21*= -l
            f22*= l
            Nastran_matrix_0_deg_11 = sp.expand(sp.poly(f11).coeffs()[1 if len(sp.poly(f11).coeffs()) == 2 else 0])
            Nastran_matrix_0_deg_12 = sp.expand(sp.poly(f12).coeffs()[1 if len(sp.poly(f12).coeffs()) == 2 else 0])
            Nastran_matrix_0_deg_21 = sp.expand(sp.poly(f21).coeffs()[1 if len(sp.poly(f21).coeffs()) == 2 else 0])
            Nastran_matrix_0_deg_22 = sp.expand(sp.poly(f22).coeffs()[1 if len(sp.poly(f22).coeffs()) == 2 else 0])
            Nastran_matrix_1_deg_11 = sp.expand((sp.poly(f11).coeffs()[0] if len(sp.poly(f11).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_1_deg_12 = sp.expand((sp.poly(f12).coeffs()[0] if len(sp.poly(f12).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_1_deg_21 = sp.expand((sp.poly(f21).coeffs()[0] if len(sp.poly(f21).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_1_deg_22 = sp.expand((sp.poly(f22).coeffs()[0] if len(sp.poly(f22).coeffs()) == 2 else 0))/sp.I
            Nastran_matrix_0_deg = Matrix([[Nastran_matrix_0_deg_11, Nastran_matrix_0_deg_12],[Nastran_matrix_0_deg_21, Nastran_matrix_0_deg_22]])
            Nastran_matrix_1_deg = Matrix([[Nastran_matrix_1_deg_11, Nastran_matrix_1_deg_12],[Nastran_matrix_1_deg_21, Nastran_matrix_1_deg_22]])
    ##        print(f"{L_z = }; {L_z_d = }; {L_t = }; {L_t_d = }; {M_z = }; {M_z_d = }; {M_t = }; {M_t_d = }; ")
            Force_matrix = Matrix([[f11,f12],[f21,f22]])
    ##        print(f"{Force_matrix=}")
            Final_matrix = -w**2*M_matrix+Kk_matrix+Force_matrix
            eq_det = sp.Eq(Final_matrix.det(),0)
            res_det = Solve_nastran(Nastran_matrix_1_deg,Nastran_matrix_0_deg)
            res_det = [i/sp.I for i in res_det]
    ##        print(f"\n{eq_det=}")
    ##        print(f"\n{res_det=}\n\n")
##            print(f"{w_j=}")
            if ind>1 and abs(w_j-w_j_prev)<=0.001:
##                print(f" ")
                w_j_final.append(w_j)
                if im(w_j)<0:
                    flutter_flag = 1
                break
            if ind == 20:
                flutter_flag = 1
                w_j_final.append(w_j)
                print(f"Итерационный процесс не сходится. Принято, что взята скорость выше критической. Результат двух последних итераций:\n{w_j_prev}\n{w_j}  (Флаттер - отриц. мнимая часть)")
                break
    if len(w_j_final)<2:
        print("Схождение к одному корню")
    return flutter_flag

def get_crit_speed_both(c1, m, h, b, o, dcy, ro, l, a, Sound_speed):
    create_bdf_file(c1)
    crit_speed_Theodorsen = 0
    crit_speed_Possio = 0
    v = 2**10
    for i in range(9, -1,-1):
        print(f"{i}...      {v = }")
        if get_Theodorsen(v, c1, m, h, b, o, dcy, ro, l, a, Sound_speed) == 1:
            v -= 2**i
        else:
            v += 2**i
    if get_Theodorsen(v, c1, m, h, b, o, dcy, ro, l, a, Sound_speed) == 1:
        crit_speed_Theodorsen = v
    else:
        crit_speed_Theodorsen = v+1
    print("Проверка результата")
    if get_Theodorsen(crit_speed_Theodorsen, c1, m, h, b, o, dcy, ro, l, a, Sound_speed) == get_Theodorsen(crit_speed_Theodorsen-1, c1, m, h, b, o, dcy, ro, l, a, Sound_speed):
        print("Ошибка, найдена неверная скорость по Теодорсену")
    else:
        print(f"Критическая скорость по Теодорсену равна {crit_speed_Theodorsen}")
        
    v = 2**10
    for i in range(9, -1,-1):
        print(f"{i}...      {v = }")
        if get_Possio(v, c1, m, h, b, o, dcy, ro, l, a, Sound_speed) == 1:
            v -= 2**i
        else:
            v += 2**i
    if get_Possio(v, c1, m, h, b, o, dcy, ro, l, a, Sound_speed) == 1:
        crit_speed_Possio = v
    else:
        crit_speed_Possio = v+1
    print("Проверка результата")
    if get_Possio(crit_speed_Possio, c1, m, h, b, o, dcy, ro, l, a, Sound_speed) == get_Possio(crit_speed_Possio-1, c1, m, h, b, o, dcy, ro, l, a, Sound_speed):
        print("Ошибка, найдена неверная скорость по Поссио")
    else:
        print(f"Критическая скорость по Поссио равна {crit_speed_Possio}")
    return (crit_speed_Theodorsen, crit_speed_Possio)

        
if __name__ == "__main__":
    c1=0.09; m=2*(10**-8); h=10; b=100; o=-20; dcy=6.283; ro=1.225e-12; l=4000; a=25
    Sound_speed = 340.294*1000
    file = open("DMIG_file.pch", 'w')
    file.close()
    get_Theodorsen(300, c1, m, h, b, o, dcy, ro, l, a, Sound_speed)
##    crit_speed_both1 = get_crit_speed_both(0.09, m, h, b, o, dcy, ro, l, a, Sound_speed)
##    crit_speed_both2 = get_crit_speed_both(0.5, m, h, b, o, dcy, ro, l, a, Sound_speed)
##    crit_speed_both3 = get_crit_speed_both(1, m, h, b, o, dcy, ro, l, a, Sound_speed)
##    print(crit_speed_both1)
##    print(crit_speed_both2)
##    print(crit_speed_both3)

##    for i in range(308,322):
##        get_Theodorsen(i, c1, k2, m, h, b, o, dcy, ro, l, a, Sound_speed)
##        print(get_Possio(i, c1, k2, m, h, b, o, dcy, ro, l, a, Sound_speed))
    
##        flutter_speed = theor_calc(i, flutter_speed)
##        print(f"\n\nspeed = {i}, flutter = {flutter_speed}\n")
