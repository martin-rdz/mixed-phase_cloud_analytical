#!/usr/bin/env python

"""
implement the variables, that were originally used by Korolev-Mazin 2003 to simplify the analytical treatment of mixed-phase clouds

Author: radenz@tropos.de
"""


import numpy as np
from types import SimpleNamespace

# some constants
#

cp_d  = 1005
cp_v  = 1850 
Rg   = 8.317 #korolev
Rg = 8.31446261815324 # J K-1 mol-1 #better
# Mma  = 28.96e-3
# Mma  = 0.02896439 # that one produces 287.058, composition from the 70s
Mma  = 0.028966 # for the new century [Gatley, Herrmann, Kretzschmar 2008]
Mmv  = 0.01806  # kg mol-1
Ra   = Rg/Mma
Rv   = Rg/Mmv

g    = 9.81
rho_i= 920
rho_w= 1000.0

T0 = 273.15
Tst = 373.15

p0 = 1e5
# p0 = 1013.25e2


def density_supercooled_water(T):
    """D.  E.  Hare and C.  M.  Sorensen: Density of supercooled water 1987 JChemPhys
    valid between -33.4 to -5/+10
	
	Args:
		T: in [K]
    """
    
    a = [0.99986, 6.690e-5, -8.486e-6, 1.518e-7, -6.9484e-9, -3.6449e-10, -7.497e-12]
    t = T-T0
    rho = a[0] + a[1]*t + a[2]*t**2 + a[3]*t**3 + a[4]*t**4 + a[5]*t**5 + a[6]*t**6
    return rho*1000


def Ew_LoweFicke(T):
    """saturation vapor pressure over water, Lowe and Ficke 1974 (PK A4.1)
	
	Args:
		T: in [K]
	"""
    t = T-T0
    w0=6.107799961; w1=4.436518521e-1; w2=1.428945805e-2; w3=2.650648471e-4; 
    w4=3.031240396e-6; w5=2.034080948e-8;  w6=6.136820929e-11;
    Ew = (w0+t*(w1+t*(w2+t*(w3+t*(w4+t*(w5+t*w6))))))*100
    return Ew

def Ei_LoweFicke(T):
    """saturation vapor pressure over ice, Lowe and Ficke 1974 (PK A4.1)
	
	Args:
		T: in [K]
	"""
    t = T-T0
    i0=6.109177956; i1=.503469897; i2=1.886013408e-2; i3=4.176223716e-4; 
    i4=5.824720280e-6; i5=4.838803174e-8; i6=1.838826904e-10;  
    Ei = (i0+t*(i1+t*(i2+t*(i3+t*(i4+t*(i5+t*i6))))))*100
    return Ei

#parametrizations store
pstore = {
    "Lw": {
        "k": lambda T: 2.495e6-2.3e3*(T-T0),
        # Rodgers and Yau, cited in wikipedia
        "RY": lambda T: 1000.0 * (2500.8 - 2.36*(T-T0) + 0.0016*(T-T0)**2 - 0.00006*(T-T0)**3),
        # Fleagle and Businger (1980. p. 113) less agreement with wmo table
        "FB": lambda T: (25-0.02274*(T-T0))*1e5,
    },
    "Li": {
        "k": lambda T: -2.7273*(T-T0)**2-207.2727*(T-T0)+2.8351e+6,
        # Rodgers and Yau, cited in wikipedia
        "RY": lambda T: 1000.0 * (2834.1 - 0.29*(T-T0) - 0.004*(T-T0)**2),
        # Fleagle and Businger (1980. p. 113) less agreement with wmo table
        "FB": lambda T: (28.34 - 0.00149*(T-T0))*1e5
    },
    "D": { #water vapor diffusion in air
        # korolev default
        "k": lambda *x: (2.26e-5+1.5e-7*(x[0]-T0))*(p0/x[1]),
        # unknown default
        "?": lambda *x: 21.2e-6 * (1.0 + 0.0071*(x[0]-T0)),
        # Hall-Pruppacher 1976 (J.Atmos.Sci.)
        "HP": lambda *x: 0.211*(x[0]/T0)**1.94*(1013.25e2/x[1])*1e-4, # converting from cm2 s-1 to m2 s-1
        # VDI
        "VDI": lambda *x: 2.252/(p0)*(x[0]/T0)**1.81,
        "fixed": lambda *x: 1.69e-5
    },
    "k": { # heat (or thermal) conductivitiy in air
        "k": lambda *x: 2.424e-2+7.95e-5*(x[0]-T0),
        # unknown default
        "?": lambda *x: (1.5207e-11 * x[0]**3 ) - ( 4.8574e-8 * x[0]**2 ) + ( 1.0184e-4 * x[0] ) - ( 3.9333e-4 ),
        # Hilsenrath 1960
        "H": lambda *x: 2.411e-2 * (1+3.309e-3*(x[0]-T0) - 1.441e-6*(x[0]-T0)**2),
        # PK 1997
        "PK": lambda *x: (5.69 + 0.017*(x[0]-T0))*1e-5*1e2*4.18684,
    },
    'cp_d': {
        "fixed": lambda *x: cp_d,
        # Hilsenrath 1960 (Tables of Thermodynamic Transport Properties..., Pergamon Press)
        "H": lambda *x: 1005.60 + 0.017211*(x[0]-T0)+0.000392*(x[0]-T0)**2
    },
    'cp_v': {
        "fixed": lambda *x: cp_v,
        # Reid 1987 (The Properties of Gases and Liquids, McGraw-Hill)
        "R": lambda *x: 1858 + 3.820e-1*(x[0]-T0)+4.220e-4*(x[0]-T0)**2-1.996e-7*(x[0]-T0)**3
    },
    "Ew": { #saturation pressure water
        # Magnus
        "LF": lambda T: Ew_LoweFicke(T),
        "M": lambda T: 100 * 6.1078 * np.exp(17.269388  * (T-T0) / (T - 35.86)),
        # Arden Buck equation
        "AB": lambda T: 6.1121e2 * np.exp((18.678-(T-T0)/234.5)*((T-T0)/(257.14+T-T0))),
        # Goff–Gratch equation
        "GG": lambda T: 10**(-7.90298*(Tst/T-1)+5.02808*np.log10(Tst/T) \
              - 1.3816e-7*(10**(11.344*(1-(T/Tst)))) \
              + 8.1328e-3 * (10**(-3.49149*(Tst/T - 1))-1) \
              + np.log10(1013.25e2))
    },
    "Ei": { #saturation pressure water
        # Magnus
        "LF": lambda T: Ei_LoweFicke(T),
        "M": lambda T: 100 * 6.1078 * np.exp(21.8745584 * (T-T0) / (T - 7.66 )),
        # Arden Buck equation
        "AB": lambda T: 6.1115e2 * np.exp((23.036-(T-T0)/333.7)*((T-T0)/(279.82+T-T0))),
        # Goff–Gratch equation
        "GG": lambda T: 10**(-9.09718*(T0/T-1) - 3.56654*np.log10(T0/T) \
              +0.876793*(1-T/T0)+np.log10(6.1071e2))
    },
    "rho_w": { # density of water
        # unknown default
        "fixed": lambda *x: rho_w,
        #
        "HS": density_supercooled_water,
    },
	"eta_a": { # dynamic viscosity of air
        # Hilsenrath 1960
        "H": lambda *x: (1.458e-6*x[0]**(3/2))/(x[0]+110.4),
        #
        "fixed": lambda *x: 1.8e-5, #kg/ ms,
    },
}

def first_key(d):
    """get the key of the first entry as a default"""
    return list(d.keys())[0]


def get_factors(
    T, p, rh='ice', c=1, verbose=False, **kwargs):
    """get the factors a0, a1, a2, Aw, Ai, ... as used in the papers
    
    Args:
        T: temperature in [K]
        p: pressure in [hPa]
        rh: either `'ice'` or relative humidity over water in [%]
        c: ice growth capacitance
        verbose: (default `False`) some diagnostic output
        **kwargs: select different parametrizations
		
    Returns:
        a SimpleNamespace with the factors
    """
    
    print('input kwargs: ', kwargs)
    f = SimpleNamespace()

    f.T = T
    f.p = p
    f.Lw_param = kwargs.get('Lw', first_key(pstore['Lw']))
    f.Lw = pstore['Lw'][f.Lw_param](T)
    f.Li_param = kwargs.get('Li', first_key(pstore['Li']))
    f.Li = pstore['Li'][f.Li_param](T) 
    
    f.D_param = kwargs.get('D', first_key(pstore['D']))
    f.D = pstore['D'][f.D_param](T, p)
    f.k_param = kwargs.get('k', first_key(pstore['k']))
    f.k = pstore['k'][f.k_param](T)
    
    f.Ew_param = kwargs.get('Ew', first_key(pstore['Ew']))
    f.Ew = pstore['Ew'][f.Ew_param](T)
    f.Ei_param = kwargs.get('Ei', first_key(pstore['Ei']))
    f.Ei = pstore['Ei'][f.Ei_param](T)
    
    f.rho_w_param = kwargs.get('rho_w', first_key(pstore['rho_w']))
    f.rho_w = pstore['rho_w'][f.rho_w_param](T)
	
    f.eta_a_param = kwargs.get('eta_a', first_key(pstore['eta_a']))
    f.eta_a = pstore['eta_a'][f.eta_a_param](T)
    
    # initial supersaturation over ice
    if rh == 'ice':
        sup0 = 0 # fraction not percent
        f.e0 = f.Ei*(1+sup0)
    else:
        f.e0 = rh*f.Ew/100.
    
    ma = 1
    mv=ma*Mmv*f.e0/(Mma*(p-f.e0))
    print('mv ', mv) if verbose else None
    mt = ma + mv
    print('weigth of adiabatic blob ', mt) if verbose else None
    f.Rt = Rg*(mv*Mma+ma*Mmv)/((ma+mv)*Mma*Mmv)
    print('gas constant moist air ', f.Rt) if verbose else None
    
    f.cp_d_param = kwargs.get('cp_d', first_key(pstore['cp_d']))
    cp_d = pstore['cp_d'][f.cp_d_param](T)
    f.cp_v_param = kwargs.get('cp_v', first_key(pstore['cp_v']))
    cp_v = pstore['cp_v'][f.cp_v_param](T)
    
    f.cp = (ma*cp_d+mv*cp_v)/(ma+mv)
    f.mv = mv
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # somehow fixed value in dQ_vs_Nr_P_T_mixed.m
    #f.cp = 1006.1
    
    qv_abd = (f.e0/(Rv*T))/((p-f.e0)/(Ra*T))
    print("qv", qv_abd, mv) if verbose else None
    print('cpt vs cp', f.cp, (1 + 0.87*qv_abd)*1003) if verbose else None
    
#     print('before', qv, T)
#     T += T*0.61*qv
#     qv   = (f.Ew/(Rv*T))/((p-f.Ew)/(Ra*T))
#     print('after', qv, T)
    
    f.rho_a= p/(f.Rt*T)
    f.rho_dry = p/(Ra*T)
    print('rho_air (with Rt)', p/(f.Rt*T), '\n Ra', p/(Ra*T)) if verbose else None
    f.ksi = f.Ew/f.Ei
	
    f.nu_a = f.eta_a/f.rho_dry

    # correction factors for the pinsky notation
    f.upsilon_w =((3*f.rho_a)/(4*np.pi*f.rho_w))**(1/3)
    f.upsilon_i = ((3*f.rho_a)/(4*np.pi*rho_i))**(1/3)/(f.ksi*c)
    
    f.a0 = g/(f.Rt*T)*(f.Lw*f.Rt/(f.cp*Rv*T)-1)
    f.a1 = 1.0/mv + (f.Lw*f.Lw)/(f.cp*Rv*T**2)
    f.a2 = 1.0/mv + (f.Lw*f.Li)/(f.cp*Rv*T**2)
    f.a3 = 1.0/mv + (f.Li*f.Li)/(f.cp*Rv*T**2)
	
	# Growth factors for ice and water
	#
    # exactly what is written in the Korolev papers
    f.Aw_alt = 1.0/((f.rho_w*f.Lw*f.Lw/(f.k*Rv*T**2))+(f.rho_w*Rv*T/(f.Ew*f.D)))
    f.Ai_alt = 1.0/((rho_i*f.Li*f.Li/(f.k*Rv*T**2))+(rho_i*Rv*T/(f.Ei*f.D)))
    # different coefficients from the matlab implementation
	# the -1 is frequently omitted (Lamb Verlinde Eq. 8.18, Khvorostyanov Curry Eq. 5.2.15b)
	# ... slight differences at all temperatures
    rho_v = f.e0/(Rv*T)
    kw2 = (f.Lw/(Rv*T)-1)*f.Lw*rho_v*f.D/(T*f.k)
    kw1 = f.D*rho_v/(f.rho_w*(kw2+1))
    f.Aw = kw1
    ki2 = (f.Li/(Rv*T)-1)*f.Li*rho_v*f.D/(T*f.k)
    ki1 = c*f.D*rho_v/(rho_i*(ki2+1))
    f.Ai = ki1
    
    f.Bi_s = 4*np.pi*rho_i*(f.ksi-1)*c*f.Ai/f.rho_a
    f.Bi   = 4*np.pi*rho_i*f.ksi*c*f.Ai/f.rho_a
    f.Bi_0 = 4*np.pi*rho_i*c*f.Ai/f.rho_a
    
    f.bi_s = f.a2*f.Bi_s
    f.bi_0 = f.a3*f.Bi_0
    f.eta  = f.a2*f.Bi_0/f.a0
    
    f.Bw   = 4*np.pi*f.rho_w*f.Aw/f.rho_a
    f.bw   = f.a1*f.Bw
    f.bi   = f.a2*f.Bi
    
    return f
    
    
# print('-- -15 -------------------------')
# factors = get_factors(T0-15, 800e2, verbose=True)  

# print('-- -35 -------------------------')
# factors = get_factors(T0-35, 800e2, verbose=True)  

# print(factors.__dict__.keys())