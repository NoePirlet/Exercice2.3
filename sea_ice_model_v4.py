import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


#definition des constantes :

#Physical constants
Lfus = 3.35e5        #Latent heat of fusion for water [J/kg]
rhoi = 917           #Sea ice density [kg/m3]
ki = 2.2             #Sea ice thermal conductivity [W/m/K]
ks = 0.31            #Snow thermal conductivity [W/m/K]
sec_per_day = 86400  #Second in one day [s/day]

#Bottom boundary condition
T_bo = -1.8          #Bottom temperature [C]



#conditions initiales
thick_0 = 0.1      #épaisseur initiale [m]
T_0 = -10           #température initiale [C]
Q_w = 2            #flux de chaleur océanique [W/m2]
snow_0 = 0.00       #initiale snow thickness [m]
nbry = 1           #nombre d'année souhaité


N_d = nbry*365  #nombre de jour dans x années

doy = (np.arange(0,N_d))%365+1 #je mets pas le +1 car je pense que c'est pour matlab ? sisi il faut le mettre sinon 364j dans l'année :/ enfaite si je le met pas c'est comme si le premier jour c'est le jour 0

#Surface heat fluxes-----------------------------------------------------------------------------------
#Other physical constants
epsilon = 0.99          #surface emissivity
sigma = 5.67e-8         #Stefan-Boltzmann constantes
Kelvin = 273.15         #Conversion from Celsius to Kelvin
alb = 0.8              #surface albedo


Q_SOL = np.zeros(N_d)       #flux de chaleur solaire
Q_NSOL = np.zeros(N_d)      #flux de chaleur non solaire

def solar_flux(day):
    
    Q_sol= 314 * np.exp((-(day-164)**2)/4608)
    return(Q_sol)

def non_solar_flux(day):
    
    Q_nsol = 118 * np.exp((-0.5*(day-206)**2)/(53**2))+179
    return(Q_nsol)
    
    
#Surface temperature----------------------------------------------------------------------------------
T_0K = T_0 + Kelvin        #convertion en Kelvin
T_boK = T_bo + Kelvin      #convertion en Kelvin
x0 = 200.15                #valeur de départ N-R

def f(T_su,h,ki,epsilon,sigma,alb,day):
    return (-(h*epsilon*sigma)/ki*T_su**4 - T_su + h/ki*((solar_flux(day))*(1-alb)+non_solar_flux(day))+T_boK)
 
def df(T_su,h,ki,epsilon,sigma,alb,day):
    return (-(4*h*epsilon*sigma)/ki*T_su**3-1)


def get_root(h,ki,epsilon,sigma,alb,day):
    root = optimize.newton(f,x0,fprime = df, args= (h,ki,epsilon,sigma,alb,day))
    return(root)
    
    
#Couple temperature and thickness---------------------------------------------------------------------


thick_temp = np.zeros(N_d) 
thick_temp[0] = thick_0
Q_csurf = np.zeros(N_d)
ROOT_surf = np.zeros(N_d)
Q_net = np.zeros(N_d)


def get_thick_temp(T_boK, ki, ks, rhoi, Lfus, snow_0, Q_w, epsilon,sigma, alb, N_d):

    for j in range(0,N_d-1):
        doy = (np.arange(0,N_d))%365+1
        ROOT_surf[j] = get_root(thick_temp[j],ki,epsilon,sigma,alb,doy[j])
        
        
        if ROOT_surf[j] >= 273.15:
            ROOT_surf[j] = 273.15
            Q_net[j] = (1-alb)*solar_flux(doy[j]) + non_solar_flux(doy[j]) - epsilon*sigma*ROOT_surf[j]**4 
            Q_csurf[j] = (ki*ks*(thick_temp[j]+snow_0))/(ki*snow_0+ks*thick_temp[j]) * (ROOT_surf[j]-T_boK)/(thick_temp[j]+snow_0)    
            thick_temp[j+1] = thick_temp[j] - (Q_csurf[j]+Q_w+Q_net[j])*sec_per_day/(rhoi*Lfus)
            
        else:    
            Q_csurf[j] = (ki*ks*(thick_temp[j]+snow_0))/(ki*snow_0+ks*thick_temp[j]) * (ROOT_surf[j]-T_boK)/(thick_temp[j]+snow_0)    
            thick_temp[j+1] =  thick_temp[j] - (Q_csurf[j]+Q_w)*sec_per_day/(rhoi*Lfus)
        
    return(thick_temp, ROOT_surf) 


#Sortie des tableaux et figures-----------------------------------------------------------------------

thick_temp = get_thick_temp(T_boK, ki, ks, rhoi, Lfus, snow_0, Q_w, epsilon,sigma, alb, N_d)[0]  #vecteur épaisseur de la glace en prenant en compte la variation de température en surface

ROOT_surf[-1] = get_root(thick_temp[-1], ki, epsilon,sigma,alb,doy[-1])   #attention ca semble donner une mauvaise valeur car il y a un saut entre la valeur fin-1 et fin :/
Q_csurf[-1] = (ki*ks*(thick_temp[-1]+snow_0))/(ki*snow_0+ks*thick_temp[-1]) * (ROOT_surf[-1]-T_boK)/(thick_temp[-1]+snow_0)

year = np.arange(1,N_d+1)


fig7 = plt.figure()
plt.plot(year, ROOT_surf)

plt.close

fig6 = plt.figure()
plt.plot(year, thick_temp)
plt.show()
plt.close


