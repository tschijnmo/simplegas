"""
Simple gas property calculation
==========================

This module contains some simple functions for calculating gas properties for
GCMC simulations.

"""

from __future__ import division

import math

import numpy as np


def peng_robinson(T,P,Tc,Pc,w,MW,Liquido):
    """Calculates gas property by Peng-Robinson Equation
    
    Calculates the compressibility factor,fugacity coefficient and density of a
    pure compound with the Peng Robinson equation of state (PR EOS).

    This function is adapted from the MATLAB function by Pinero,  R, Serna, JG,
    and Martin, A.

    Bibliography: 
      - ORBEY. H, SANDLER. I; Modeling Vapor-Liquid Equilibria: cubic equations of state and their mixing rules; Cambridge University Press (1998)
      - Walas,Stanley. M ; Phase Equilibria in Chemical Engineering ; Boston,
	Butterworth Publishers (1984)                 
    
    function result = peng_robinson(T,P,Tc,Pc,w,MW,Liquido)
    :param T: Temperature in K                                          
    :param P: Presure in Pa                                             
    :param Tc: critical temperature in K                               
    :param Pc: critical presure in Pa                                   
    :param w: accentic factor
    :param MW: molar weigth in kg/mol
    :param Liquido:  if Liquido is 1, then calculates liquid fugacity;  if
        Liquido is 0 then calculates vapor fugacity

    Example:
    Z, fhi, density = PengRobinson(273,2*1.013*1e5,304.21,7.382*1e6,0.225,0.044,1)

    """


    R = 8.314; # gas constant in J/(mol K)

    # Reduced variables
    Tr = T/Tc ;
    Pr = P/Pc ;

    # Parameters of the EOS for a pure component
    m = 0.37464 + 1.54226*w - 0.26992*w ** 2;
    alfa = (1 + m*(1 - math.sqrt(Tr))) ** 2;
    a = 0.45724*(R*Tc) ** 2/Pc*alfa;
    b = 0.0778*R*Tc/Pc;
    A = a*P/(R*T) ** 2;
    B = b*P/(R*T);

    # Compressibility factor
    #  Z = np.roots([-(A * B - B ** 2 - B ** 3), (A - 3 * B ** 2 - 2 * B), -(1 - B), 1]);
    Z = np.roots([1.0, -(1 - B), (A - 3 * B ** 2 - 2 * B), -(A * B - B ** 2 - B ** 3)]);

    ZR = [];
    for root in Z:
        if isinstance(root, complex):
            continue
        else:
            ZR.append(root)

    if Liquido == 1:
        Z = min(ZR);   
    else:
        Z = max(ZR);

    # Fugacity coefficient
    fhi = math.exp(Z - 1 - math.log(Z-B) - A/(2*B*math.sqrt(2))*math.log((Z+(1+math.sqrt(2))*B)/(Z+(1-math.sqrt(2))*B)));
    if isinstance(fhi, complex):
        raise ValueError(
            'No real solution for "fhi" is available in this phase'
            )
    density=P*MW/(Z*R*T);

    return Z, fhi, density
        

def chemical_pot(MW, fhi, T, P):
    """Computes the chemical potential

    This function just uses the statistical mechanics of the ideal gas to get
    the reference thermal wave length fot computing the chemical potential,
    then the fugacity coefficient can be applied to the density for correction.

    """

    toAA2 = 9.5725E27;
    h = 4.135667516E-15;
    kB = 8.6173324E-5;

    lamb = math.sqrt( (h ** 2 * toAA2 ) / (2 * math.pi * MW * 1000.0 * kB * T) );  # Thermal wave length
    rho = ((P * fhi) / (kB * T) ) * 6.241509E-12;
    return  kB * T * math.log( lamb ** 3 * rho );


def h2_mu(T, P):
    """Computes the chemical potential of the hydrogen gas

    This is a wrapper function for computing the chemical potential of the hydrogen gas.
    """

    Tc = 32.97;
    Pc = 1.293E6;
    w = -0.215;
    MW = 0.00201588;
    Liquido = 0;

    Z,fhi,density = peng_robinson(T,P,Tc,Pc,w,MW,Liquido)
    return chemical_pot(MW, fhi, T, P)




