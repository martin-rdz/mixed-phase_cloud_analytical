#!/usr/bin/env python

import numpy as np

TsIGR = np.array([-40. , -39.6, -39.2, -38.8, -38.4, -38. , -37.6, -37.2, -36.8,
       -36.4, -36. , -35.6, -35.2, -34.8, -34.4, -34. , -33.6, -33.2,
       -32.8, -32.4, -32. , -31.6, -31.2, -30.8, -30.4, -30. , -29.6,
       -29.2, -28.8, -28.4, -28. , -27.6, -27.2, -26.8, -26.4, -26. ,
       -25.6, -25.2, -24.8, -24.4, -24. , -23.6, -23.2, -22.8, -22.4,
       -22. , -21.6, -21.2, -20.8, -20.4, -20. , -19.6, -19.2, -18.8,
       -18.4, -18. , -17.6, -17.2, -16.8, -16.4, -16. , -15.6, -15.2,
       -14.8, -14.4, -14. , -13.6, -13.2, -12.8, -12.4, -12. , -11.6,
       -11.2, -10.8, -10.4, -10. ,  -9.6,  -9.2,  -8.8,  -8.4,  -8. ,
        -7.6,  -7.2,  -6.8,  -6.4,  -6. ,  -5.6,  -5.2,  -4.8,  -4.4,
        -4. ,  -3.6,  -3.2,  -2.8,  -2.4,  -2. ,  -1.6,  -1.2,  -0.8,
        -0.4])

IGR_planar = np.array([0.8683977 , 0.8683977 , 0.8683977 , 0.8683977 , 0.8683977 ,
       0.86803007, 0.8675514 , 0.86707272, 0.86659405, 0.86611538,
       0.86563671, 0.86515804, 0.86467937, 0.8642007 , 0.86372203,
       0.86324336, 0.86284572, 0.86284572, 0.86284572, 0.86284572,
       0.86284572, 0.86284572, 0.86284572, 0.86284572, 0.86284572,
       0.86284572, 0.86284572, 0.86284572, 0.86284572, 0.86284572,
       0.86284572, 0.86284572, 0.86284572, 0.86284572, 0.86284572,
       0.86284572, 0.86284572, 0.86284572, 0.86284572, 0.86284572,
       0.86284572, 0.86284572, 0.86284572, 0.86284572, 0.86284572,
       0.86284572, 0.86284572, 0.86284572, 0.86284572, 0.84964364,
       0.83027136, 0.78624718, 0.74120905, 0.69764403, 0.6574929 ,
       0.6187466 , 0.5940228 , 0.57246002, 0.55511304, 0.5383176 ,
       0.52890477, 0.52174049, 0.51748596, 0.51652568, 0.52007336,
       0.52787712, 0.54545168, 0.56690622, 0.59418496, 0.62953437,
       0.67341242, 0.72370159, 0.77535018, 0.83264917, 0.88646885,
       0.94006683, 0.99712498, 1.05191196, 1.10876756, 1.18442555,
       1.25482499, 1.34347271, 1.42482849, 1.48089908, 1.505418  ,
       1.50140197, 1.44834316, 1.32614029, 1.19117235, 1.00422648,
       0.81398546, 0.70675361, 0.72078798, 0.80369398, 0.87167576,
       0.89329019, 0.91009479, 0.92555973, 0.9448907 , 0.96843934])

IGR_col = np.array([1.22392513, 1.21148437, 1.20116628, 1.19363228, 1.18433045,
       1.17232076, 1.16031107, 1.14645447, 1.13522645, 1.13009043,
       1.1249544 , 1.12074841, 1.11747831, 1.11457384, 1.11250617,
       1.1104385 , 1.10837083, 1.10912482, 1.11035424, 1.11158366,
       1.11281309, 1.11404251, 1.11520257, 1.11520257, 1.11520257,
       1.11520257, 1.11520257, 1.11427283, 1.11245328, 1.11063373,
       1.10881418, 1.10714927, 1.10559072, 1.10403218, 1.10247363,
       1.10092626, 1.09960541, 1.09828456, 1.0969637 , 1.09564285,
       1.094322  , 1.09394935, 1.09394935, 1.09394935, 1.0879664 ,
       1.07983397, 1.06885591, 1.05165664, 1.01879261, 0.92117849,
       0.83561374, 0.78624718, 0.74120905, 0.69764403, 0.6574929 ,
       0.6187466 , 0.5940228 , 0.57246002, 0.55511304, 0.5383176 ,
       0.52890477, 0.52174049, 0.51748596, 0.51652568, 0.52007336,
       0.52787712, 0.54545168, 0.56690622, 0.59418496, 0.62953437,
       0.67341242, 0.72370159, 0.77535018, 0.83264917, 0.88646885,
       0.94006683, 0.99712498, 1.05191196, 1.10876756, 1.18442555,
       1.25482499, 1.34347271, 1.42482849, 1.48089908, 1.505418  ,
       1.50140197, 1.44834316, 1.32614029, 1.19117235, 1.00422648,
       0.81398546, 0.70675361, 0.72078798, 0.80369398, 0.87167576,
       0.89329019, 0.91009479, 0.92555973, 0.9448907 , 0.96843934])

IGR_chen = np.array([1.22382494, 1.22382494, 1.22382494, 1.22382494, 1.22382494,
       1.22382494, 1.22382494, 1.22382494, 1.22382494, 1.22382494,
       1.22382494, 1.22382494, 1.22382494, 1.22382494, 1.22382494,
       1.22382494, 1.22382494, 1.22382494, 1.22382494, 1.22382494,
       1.22382494, 1.22382494, 1.22382494, 1.22382494, 1.22382494,
       1.22382494, 1.24566689, 1.26930506, 1.28759672, 1.30293771,
       1.3182787 , 1.35067257, 1.38477743, 1.42004963, 1.45774378,
       1.49583913, 1.54233309, 1.60178596, 1.67586763, 1.75307581,
       1.8331452 , 1.87195614, 1.89602349, 1.89797522, 1.81147807,
       1.64774019, 1.43345859, 1.23670029, 1.06494624, 0.91395885,
       0.78058501, 0.68395666, 0.59660315, 0.52502151, 0.47129959,
       0.42196933, 0.37632453, 0.34356317, 0.31967118, 0.30183852,
       0.28885975, 0.27872854, 0.27341152, 0.2732915 , 0.27940562,
       0.29417697, 0.3147245 , 0.34243901, 0.37818303, 0.418501  ,
       0.46944331, 0.53080276, 0.59667274, 0.67489539, 0.76759512,
       0.87138387, 0.99055847, 1.12285885, 1.27123386, 1.46897874,
       1.66012123, 1.84799253, 2.04610996, 2.21810007, 2.31239406,
       2.29497846, 2.14873942, 1.80103869, 1.39754246, 0.93880334,
       0.59530584, 0.47851276, 0.59747185, 0.74558773, 0.79310236,
       0.83609253, 0.87292277, 0.90536772, 0.93421052, 0.96305332])

IGR_H07 = np.array([1.11695309, 1.11695309, 1.11695309, 1.11695309, 1.11695309,
       1.11695309, 1.11695309, 1.11695309, 1.11695309, 1.11695309,
       1.11695309, 1.11695309, 1.11695309, 1.11695309, 1.11695309,
       1.11695309, 1.11695309, 1.11695309, 1.11695309, 1.11695309,
       1.11695309, 1.11695309, 1.11695309, 1.11695309, 1.11695309,
       1.11695309, 1.12059596, 1.12600006, 1.13372335, 1.14163724,
       1.15017637, 1.15950318, 1.16986433, 1.18022548, 1.1944249 ,
       1.21431414, 1.23420337, 1.257119  , 1.28264275, 1.30816651,
       1.33341982, 1.35841084, 1.36247018, 1.35782687, 1.32238451,
       1.27577643, 1.21055163, 1.13563121, 1.05489794, 0.97399127,
       0.89341442, 0.82459275, 0.76268059, 0.70929863, 0.66502801,
       0.62721455, 0.59754828, 0.57448397, 0.55511993, 0.54235262,
       0.53275046, 0.52357949, 0.52172164, 0.52088857, 0.52497182,
       0.53442841, 0.55037185, 0.56820695, 0.59593785, 0.63140045,
       0.67583744, 0.72562326, 0.77952022, 0.83494363, 0.89008008,
       0.94449616, 0.99837937, 1.05204312, 1.11526206, 1.18084031,
       1.26301658, 1.34623255, 1.43097584, 1.48506323, 1.51882874,
       1.51153131, 1.45381721, 1.33414712, 1.18325023, 0.99877254,
       0.81078845, 0.69863504, 0.71857491, 0.8038081 , 0.87415236,
       0.89223321, 0.91223881, 0.93266111, 0.95389585, 0.97028555])


rho_dep_Miller = np.array([0.9    , 0.9    , 0.9    , 0.9   , 0.9,
       0.9       , 0.89202734, 0.88405468, 0.87608202, 0.86810936,
       0.8601367 , 0.85216404, 0.84419138, 0.83621872, 0.82824606,
       0.8202734 , 0.81245476, 0.80657581, 0.80069687, 0.79481793,
       0.78893898, 0.78306004, 0.77399318, 0.76414963, 0.75431581,
       0.7447535 , 0.73519118, 0.72519532, 0.71511168, 0.70532358,
       0.69641751, 0.68751143, 0.67885124, 0.67034999, 0.66184874,
       0.65277261, 0.6435853 , 0.63580693, 0.62810328, 0.62044377,
       0.61355328, 0.60666279, 0.60082869, 0.59590691, 0.59127082,
       0.58689591, 0.58238073, 0.57745895, 0.57253717, 0.56784696,
       0.56390954, 0.51165512, 0.44182156, 0.37743443, 0.32291473,
       0.27708721, 0.23413436, 0.20076732, 0.17341801, 0.15165737,
       0.13197026, 0.11871902, 0.11198575, 0.11095493, 0.11703477,
       0.13119579, 0.15104038, 0.18667795, 0.23767038, 0.30755112,
       0.39439   , 0.46862608, 0.55464684, 0.63638786, 0.72018635,
       0.80368649, 0.88606636, 0.89791526, 0.89874535, 0.89950255,
       0.89908096, 0.89842472, 0.8228276 , 0.74195322, 0.66088393,
       0.58741223, 0.53277055, 0.50903478, 0.50847584, 0.52331681,
       0.64141884, 0.88713135, 0.89926831, 0.90020579, 0.90086741,
       0.90086741, 0.90086741, 0.90058086, 0.90026332, 0.89994579])



def get_igr(T, source):
    """return the IGR form the data
    
    C94 is also what should use ISMAHEL
    
    Args:
        T: temperature [K]
        source: string H08col, H08plan, C94, H07
    """
    assert np.all(T > 0)

    if source=='H08col':
        igr_here = np.interp(T-273.15, TsIGR, IGR_col)
    elif source=='H08plan':
        igr_here = np.interp(T-273.15, TsIGR, IGR_planar)
    elif source=='C94':
        igr_here = np.interp(T-273.15, TsIGR, IGR_chen)
    elif source=='H07':
        igr_here = np.interp(T-273.15, TsIGR, IGR_H07)
    else:
        raise ValueError()
    
    return igr_here

def get_rho_dep(T):
    """get the deposition density for Miller 1979
    
    
    """
    return np.interp(T-273.15, TsIGR, rho_dep_Miller)*1e3
