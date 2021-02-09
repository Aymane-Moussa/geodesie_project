import math
from mpmath import *

def calcul_de_parametre(a,invf):
    f = 1 / invf
    b = a - a * f
    e = math.sqrt((pow(a,2) - pow(b,2)) / pow(a,2))
    ep = math.sqrt((pow(a,2) - pow(b,2)) / pow(b,2))
    alpha = math.acos(1 - f)
    c = pow(a,2) / b
    return e,ep,alpha,c

def dms2degrees(L):
    if(L[0]<0):
        return  L[0] - (L[1]/60) - (L[2]/3600)
    else:
        return L[0] + (L[1]/60) + (L[2]/3600)


def decdeg2dms(dd):
   is_positive = dd >= 0
   dd = abs(dd)
   minutes,seconds = divmod(dd*3600,60)
   degrees,minutes = divmod(minutes,60)
   degrees = degrees if is_positive else -degrees
   return degrees,minutes,seconds

def calcul_rayon_courbure(a, invf, φ):
    angleInDegrees = dms2degrees(φ)
    e= calcul_de_parametre(a,invf)
    w = math.sqrt(1 - pow((e[0]*math.sin(math.radians(angleInDegrees))),2))
    M = a * (1 - pow(e[0],2)) / pow(w,3)
    N = a / w
    return N,M,w

def calcul_rayon_selon_azimut_alpha(a,invf,φ,alpha):
    M=calcul_rayon_courbure(a, invf, φ)
    N=calcul_rayon_courbure(a, invf, φ)
    R_alpha= (M[1]*N[0])/(M[1]*pow(math.sin(math.radians(alpha)),2)+N[0]*pow(math.cos(math.radians(alpha)),2))
    return R_alpha

def Transformation_coord_geo_vers_cart(φ,λ,h,invf,a):
    angleInDegrees = dms2degrees(φ)
    angleInDegrees2 = dms2degrees(λ)
    N=calcul_rayon_courbure(a, invf, φ)
    X = (N[0]+ h) * (math.cos(math.radians(angleInDegrees))) * math.cos(math.radians(angleInDegrees2))
    Y = (N[0] + h) * math.cos(math.radians(angleInDegrees)) * math.sin(math.radians(angleInDegrees2))
    e= calcul_de_parametre(a, invf)
    Z = (N[0] * (1 - pow(e[0],2)) + h) * math.sin(math.radians(angleInDegrees))
    return X,Y,Z

def Transformation_coord_cart_vers_geo_direct(X,Y,Z,invf,a):
    e = calcul_de_parametre(a, invf)
    λ=math.atan(Y/X)
    r = math.sqrt(pow(X,2) + pow(Y,2) + pow(Z,2))
    mu = math.atan((Z / math.sqrt(pow(X,2) + pow(Y,2))) * ((1 - 1 / invf) + (a * pow(e[0],2) / r)))
    φ = math.atan((Z * (1 - 1 / invf) + pow(e[0],2) * a * pow(math.sin(mu),3)) / ((1 - 1 / invf) * (math.sqrt(pow(X,2) + pow(Y,2)) - pow(e[0],2) * a * pow(math.cos(mu),3))))
    h = math.sqrt(pow(X,2) + pow(Y,2)) * math.cos(φ) + Z * math.sin(φ) - a * math.sqrt(1 - pow(e[0],2) * pow(math.sin(φ),2))
    return math.degrees(φ),math.degrees(λ),h

def Transformation_coord_cart_vers_geo(X,Y,Z,invf,a,precision):
    e = calcul_de_parametre(a, invf)
    φ0=math.atan((Z/math.sqrt(pow(X,2)+pow(Y,2)))*(1+(pow(e[0],2)/(1-pow(e[0],2)))))
    N0 = a / (math.sqrt(1 - pow(e[0],2) * pow(math.sin(φ0),2)))
    φ1=math.atan ((Z/math.sqrt(pow(X,2)+pow(Y,2))*(1+(N0*pow(e[0],2)*math.sin(φ0)/(1-pow(e[0],2))))))
    N1 = a / (math.sqrt(1 - pow(e[0],2) * pow(math.sin(φ0),2)))
    Nav= N1
    φap = φ1
    φav = φ0
    while abs(φap - φav)<precision :

            φav = φap
            φap = math.atan((Z / math.sqrt(pow(X,2) + pow(Y,2))) * (1 + (Nav * pow(e[0],2) * math.sin(φav) / (1 - pow(e[0],2)))))
            Nap = a / (math.sqrt(1 - pow(e[0],2) * pow(math.sin(φap),2)))
            Nav = Nap

    N = Nav
    φ = φav
    λ= math.atan(Y/X)
    h = (Z / math.sin(φ)) - (N * (1 - pow(e[0],2)))
    φ = decdeg2dms(math.degrees(φ))
    λ =decdeg2dms(math.degrees(λ))
    return φ[0],φ[1],φ[2],λ[0],λ[1],λ[2],h

def Calcul_des_latitudes_reduite_geocentrique (φ,a,b):
    φ=dms2degrees(φ)
    β=math.degrees(math.atan((b/a)*math.tan(math.radians(φ))))
    ψ= math.degrees(math.atan((pow(b,2)/pow(a,2))*math.tan(math.radians(φ))))
    β=decdeg2dms(β)
    ψ=decdeg2dms(ψ)
    return β[0],β[1],β[2],ψ[0],ψ[1],ψ[2]
def Calcul_des_latitudes_reduite_geodesique(ψ,a,b):
    ψ=dms2degrees(ψ)
    φ=math.degrees(math.atan((pow(a,2)/pow(b,2))*math.tan(math.radians(ψ))))
    β=math.degrees(math.atan((a/b)*math.tan(math.radians(ψ))))
    φ = decdeg2dms(φ)
    β = decdeg2dms(β)
    return φ[0],φ[1],φ[2],β[0],β[1],β[2]
def Calcul_des_latitudes_geodesique_geocentrique(β,a,b):
    β=dms2degrees(β)
    φ=math.degrees(math.atan((a/b)*math.tan(math.radians(β))))
    ψ=math.degrees(math.atan((b/a)*math.tan(math.radians(β))))
    φ =decdeg2dms(φ)
    ψ= decdeg2dms(ψ)
    return φ[0],φ[1],φ[2],ψ[0],ψ[1],ψ[2]

def Calcul_de_longueur_d_un_arc_de_meridien(φ,a,invf):
    e=calcul_de_parametre(a,invf)
    A=1+(3/4)*pow(e[0],2)+(45/64)*pow(e[0],4)+(175/256)*pow(e[0],6)+(11025/16348)*pow(e[0],8)+(43659/65536)*pow(e[0],10)
    B=(3/4)*pow(e[0],2)+(15/16)*pow(e[0],4)+(525/512)*pow(e[0],6)+(2205/2048)*pow(e[0],8)+(72765/65536)*pow(e[0],10)
    C=(15/16)*pow(e[0],4)+(105/256)*pow(e[0],6)+(2205/4096)*pow(e[0],8)+(10395/16348)*pow(e[0],10)
    D=(35/512)*pow(e[0],6)+(315/2048)*pow(e[0],8)+(2485/131072)*pow(e[0],10)
    E=(315/16384)*pow(e[0],8)+(3465/85536)*pow(e[0],10)
    F=(639/131072)*pow(e[0],10)
    S0_φ=a*(1-pow(e[0],2))*(A*φ-(B/2)*math.sin(2*φ)+(C/4)*math.sin(4*φ)-(D/6)*math.sin(6*φ)+(E/8)*math.sin(8*φ)-(F/10)*math.sin(10*φ))
    return S0_φ
def Calcul_de_longueur_d_un_arc_de_parallele(a,invf,φ,delta_lambda):
    e = calcul_de_parametre(a, invf)
    w = math.sqrt(1 - pow((e[0] * math.sin(φ)), 2))
    N = a / w
    Rp=N*math.cos(φ)
    L=Rp*delta_lambda
    return L
def Calcul_de_la_surface_d_une_partie_terrestre(a,invf,b,φ1,φ2,lambda1,lambda2):
    e = calcul_de_parametre(a, invf)
    Z=((lambda2-lambda1)/2)*pow(b,2)*(((1/(2*e[0]))*math.log((1+e[0]*math.sin(φ2))/(1-e[0]*math.sin(φ2)))+(math.sin(φ2))/(1-pow(e[0],2)*pow(math.sin(φ2),2))-((1/(2*e[0]))*math.log((1+e[0]*math.sin(φ1))/(1-e[0]*math.sin(φ1)))+(math.sin(φ1))/(1-pow(e[0],2)*pow(math.sin(φ1),2)))))
    return Z
def Probleme_direct(φ1,lambda1,A12,σ12):
    φ1=dms2degrees(φ1)
    lambda1=dms2degrees(lambda1)
    φ2=math.degrees(math.asin(math.sin(math.radians(φ1))*math.cos(σ12)+math.cos(math.radians(φ1))*math.sin(σ12)*math.cos(math.radians(A12))))
    lambda2=lambda1+math.degrees(acot((1/math.sin(math.radians(A12)))*(cot(σ12)*math.cos(math.radians(φ1))-math.sin(math.radians(φ1))*math.cos(math.radians(A12)))))
    A21=math.degrees(acot((1/math.sin(math.radians(A12)))*(math.cos(σ12)*math.cos(math.radians(A12))-math.tan(math.radians(φ1))*math.sin(σ12))))
    return φ2,lambda2,A21

def Proble_inverse(φ1,lambda1,φ2,lambda2):
    φ1 = dms2degrees(φ1)
    lambda1 = dms2degrees(lambda1)
    φ2 = dms2degrees(φ2)
    lambda2 = dms2degrees(lambda2)
    σ12=math.degrees(acos(math.sin(math.radians(φ1))*math.sin(math.radians(φ2))+math.cos(math.radians(φ1))*math.cos(math.radians(φ2))*math.cos(math.radians(lambda2-lambda1))))
    σ12=math.radians(σ12)
    A12=math.degrees(acot((math.tan(math.radians(φ2))*math.cos(math.radians(φ1)))/math.sin(math.radians(lambda2-lambda1))-math.sin(math.radians(φ1))*cot(math.radians(lambda2-lambda1))))
    A21=math.degrees(acot(-((math.tan(math.radians(φ1))*math.cos(math.radians(φ2)))/math.sin(math.radians(lambda2-lambda1))-math.sin(math.radians(φ2))*cot(math.radians(lambda2-lambda1)))))
    return σ12,A12,A21


#print (Transformation_coord_geo_vers_cart([33,26,59.672],[-7,33,27.295],243.730,293.46602,6378249.20))
#print(Transformation_coord_cart_vers_geo_direct(5281240.340345062,-700688.8087494428,3495569.995750439,293.46602,6378249.20))
#print(Transformation_coord_cart_vers_geo(5281240.340345062,-700688.8087494428,3495569.995750439,293.46602,6378249.20,0.01))
#print(Calcul_des_latitudes_reduite_geocentrique([33,26,59.672],6378249.20,6356515.00))
#print(calcul_rayon_courbure(6378249.2,293.730,[33,26,59.672]))
#print(calcul_rayon_selon_azimut_alpha(6378249.2,293.730,[33,26,59.672],45))
#print(Calcul_de_longueur_d_un_arc_de_meridien(math.pi/4,6378249.20,293.46602))
#print(Calcul_de_longueur_d_un_arc_de_parallele(6378249.20,293.46602,0,math.pi*2))
#print(Calcul_de_la_surface_d_une_partie_terrestre(6378249.20,293.46602,6356515.00,-1.57,1.57,0,6.28))
#print(Probleme_direct([33,26,59.672],[-7,33,27.295],22.674575,0.297785))
#print(Proble_inverse([33,26,59.672],[-7,33,27.295],[48,50,11.2],[2,20,13.8]))