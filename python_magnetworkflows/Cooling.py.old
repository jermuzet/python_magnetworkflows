#!/usr/bin/env python
# encoding: UTF-8

import math
import yaml
import json
import deserialize

# from python_magnetgeo import deserialize
from python_magnetgeo import Helix
from python_magnetgeo import InnerCurrentLead
from python_magnetgeo import OuterCurrentLead
from python_magnetgeo import Insert
import freesteam

        
##############################################################
#
#
##############################################################

def Vpump(self, current, Imax, Vpumpmax):
    """
    computes RPM for Pump as function of current
    """

    Vpumpmin = 1000
    if current > Imax:
        Rpm = Vpumpmax
    else:
        Rpm = Vpumpmin+(Vpumpmax-Vpumpmin)*(current/Imax)**2

    return Rpm

def Flow(self, Vpump, Vpumpmax, Flowmax):
    """
    computes Flow as a fonction of Vpump
    """

    Flow = Flowmax * (Vpump/Vpumpmax)
    return Flow

def PressureHP(self, BP, HP, Vpump, Vpmax):
    """
    computes HP pressure as a fonction of Vpump
    """
    return BP + (HP-BP) * (Vpump/Vpmax)**2
        
def hydraulic_diameter(self, n, Assembly):
    """
    computes the hydraulic_diameter and length of cooling channel n 
    """
        
    n0=n
    n1=1
    if n==0:
        n1=0
        n0=0
    elif n==len(Assembly.Helices):
        n1=n
        n0=0
    else:
        n1=n
        n0=n-1

    Zmin = None
    Zmax = None
    L = None
    dh = None
    S = None
    Zmin = None
    Zmax = None
        
    diameter_ratio = None
    if n==0:
        H0 = None
        with open(Assembly.Helices[n0]+".yaml", 'r') as f:
            H0 = yaml.load(f, Loader=yaml.FullLoader)

        Zmin = H0.z[0]
        Zmax = H0.z[1]

        L = abs(H0.z[1]-H0.z[0])
        e = H0.r[0]-Assembly.innerbore
        dh = 2 * e
        S = math.pi * e * ( H0.r[0]+Assembly.innerbore )
        diameter_ratio = 1 - e / Assembly.innerbore
        #print "n==0: e=", e
    elif n==len(Assembly.Helices):
        H1 = None
        with open(Assembly.Helices[n-1]+".yaml", 'r') as f:
            H1 = yaml.load(f, Loader=yaml.FullLoader)

        Zmin = H1.z[0]
        Zmax = H1.z[1]

        L = abs(H1.z[1]-H1.z[0])
        e= Assembly.outerbore-H1.r[1]    
        dh = 2 * e
        S = math.pi * e * ( Assembly.outerbore+H1.r[1] )
        diameter_ratio = 1 - e / H1.r[1]
        #print "n==%d: e="%len(self.Helices), e
    else:
        H0 = None
        with open(Assembly.Helices[n-1]+".yaml", 'r') as f:
            H0 = yaml.load(f, Loader=yaml.FullLoader)
            
        H1 = None
        with open(Assembly.Helices[n]+".yaml", 'r') as f:
            H1 = yaml.load(f, Loader=yaml.FullLoader)

        Zmin = min(H0.z[0], H1.z[0]) # Zmin of corresonding Ring
        Zmax = max(H0.z[1], H1.z[1]) # Zmax of corresponding Ring
            
        L = max(abs(H0.z[1]-H0.z[0]),abs(H1.z[1]-H1.z[0]))
        e = H1.r[0]-H0.r[1]
        dh = 2 * e
        S = math.pi * e * ( H1.r[0]+H0.r[1] )
        diameter_ratio = 1 - e / H0.r[1]
        #print "n!=0 and n!=%d: e="%len(self.Helices), e, H1.r[1], H0.r[0]
    #print "L=%g mm, a=%g" % (L,diameter_ratio)

    return (dh, L, S, Zmin, Zmax)
        
def cooling(self, Assembly, pressure=15, dpressure=25, temperature=30, correlation="colburn", friction="blasius", Pextra=1., fguess=0.055, uguess=0):
    """
    compute heat exchange coefficient

    input:
    pressure : mean pressure in bar
    delta pressure :
    temperature : mean temperature in Celsius
    correlation
    friction
    Pextra
    fguess:
    uguess:

    return:
    a tuple containing (h, Tw, flow)
    with h and Tw are list of heat exchange coeff and water temperature by channel
    """

    print ("cooling: ")
    print ("pressure=", pressure)
    print ("dpressure=", dpressure)
    print ("temperature=", temperature)
    print ("correl=", correlation)
    print ("friction=", friction)
    print ("Pextra=", Pextra)
    print ("fguess=", fguess)
    print ("uguess=", uguess)

    heat_coeff = []
    Tw = []
    flow = []
    for i in range(len(Assembly.Helices)+1):
        (dh, L, S, u, Reynolds, Prandtl, Nusselt, h, f) = self.channel(Assembly, i, pressure, dpressure, temperature, correl=correlation, friction=friction, Pextra=Pextra, fguess=fguess, uguess=uguess)
        heat_coeff.append(h)
        Tw.append(temperature)
        flow.append([dh, L, S, u, Reynolds, Prandtl, Nusselt, f])

    return (heat_coeff, Tw, flow)

def reynolds(self, dh, L, pressure, dpressure, temperature, friction="blasius", Pextra=1., fguess=0.055, uguess=0):
        
    #print "Steam(P=%g bar,T=%g C) : "%(pressure,temp)
    Steam = freesteam.steam_pT(pressure*1.e+5,temperature+273)
    rho = Steam.rho
    Cp = Steam.cp
    mu = Steam.mu
    k = Steam.k
    # print "rho=%g, Cp=%g, mu=%g, k=%g"%(rho, Cp, mu, k)

    if uguess != 0:
        Reynolds = rho*uguess*(dh*1.e-3)/mu
        f = fguess
        return (uguess, Reynolds, f)
        
    #Darcy friction factor f : Blasius, ...
    iterate=False
    f = fguess
    u = None
    while True:
        u = math.sqrt(2*(dpressure*1.e+5)/(rho*(Pextra+f*L/dh)))

        ## eventually loop over Reynolds to get f
        Reynolds = rho*u*(dh*1.e-3)/mu
        if friction == "constant":
            iterate = False
        elif friction == "blasius":
            iterate = True
            Cf = 0.3164*math.pow(Reynolds,-0.25)
        elif friction == "filonenko":
            iterate = True
            Cf = math.pow(1.82*math.log10(Reynolds)-1.64,-2)
        elif friction == "karman":
            iterate = True
            Cf  = math.pow(1.93*math.log10(Reynolds*math.sqrt(f))-0.537,-2)
        elif friction == "rough":
            iterate = True
            eps = 2.5e-2 # mm
            rstar = 1/math.sqrt(8.) * (Reynolds*math.sqrt(f))*eps/dh
            brstar = 1/(1.930*math.sqrt(f)) + math.log10(1.9/math.sqrt(8.) * eps/dh)
            ###print "brstar=%g" % brstar

            # Cf = math.pow(-1.930*math.log(1.90/(Reynolds*math.sqrt(f))*(1+0.34*rstar*math.exp(-11./rstar))),-2.)
            Cf = math.pow(-2.00*math.log10(2.51/(Reynolds*math.sqrt(f))*(1+rstar/3.3)),-2.)
            
        # Gnielinski breaks when a tends to 1
        # elif friction == "gnielinski":
        #     a = diameter_ratio
        #     Re = Reynolds * ( (1.+a**2) * math.log(a)+(1-a**2) / ( (1.-a)**2 * math.log(a) )) 
        #     Cf = math.pow(1.8*math.log10(Re)-1.5,-2)
        # # print ("%s Cf=%g" % (friction,Cf) )

        if iterate and math.fabs( (Cf-f)/f ) >= 1.e-3:
            f = Cf
        else:
            break
    return (u, Reynolds, f)
    
def channel(self, Assembly, n, pressure, dpressure, temperature, correl="colburn", friction="blasius", Pextra=1., fguess=0.055, uguess=0):
    correlation = {
        "dittus": (0.023, 0.8, 0.4),
        "colburn": (0.023, 0.8, 0.3),
        "silberberg": (0.015, 0.85, 0.3)
    }

    #friction dict to functions...
        
    if n<0 and n>len(Assembly.Helices):
        raise Exception("incoherent channel id %d"%n)

    (dh, L, S, Zmin, Zmax) = self.hydraulic_diameter(n, Assembly)
    self.Parameters["Zmin%d" % n] = Zmin * 1.e-3 # Zmin of corresonding Ring
    self.Parameters["Zmax%d" % n] = Zmax * 1.e-3 # Zmax of corresponding Ring
    self.Parameters["dh%d" % n] = dh * 1.e-3
    self.Parameters["Sh%d" % n] = S * 1.e-6

    (u, Reynolds, f) = self.reynolds(dh, L, pressure, dpressure, temperature, friction, Pextra, fguess, uguess)
        
    # Prandtl
    Steam = freesteam.steam_pT(pressure*1.e+5,temperature+273)
    Cp = Steam.cp
    mu = Steam.mu
    k = Steam.k
    Prandtl = mu * Cp / k
        
    # Nusselt= Dittus, Colburn, Silberberg,...
    Nusselt = None
    if correl in correlation.keys():
        (a, n, m) = correlation[correl]
        Nusselt = a * math.pow(Reynolds,n) * math.pow(Prandtl,m)
        h = k*Nusselt/(dh*1.e-3)
    else:
        if correl == "montgomery":
            h = 1426*(1+1.5e-2*temperature)*math.pow(u,0.8)/math.pow(dh*1.e-3,0.2)
            # print ("h=%g" % h)
        else:
            raise Exception("incoherent correlation %s\n"%correl)
        
    return (dh, L, S, u, Reynolds, Prandtl, Nusselt, h, f)

def CoolingBCs(self, Assembly, params, correlation="montgomery", friction="constant", Pextra=1., fguess=0.055):
    """
    returns a dict of Cooling Boundary conditions for given params 

    params:
    pressure: pressure on HP side
    dpressure: diffrence of pressure between HP and BP 
    temperature: input water temperature (HP side)
    uguess: water velocity
    """

    print ("CoolingsBCs: ")
    print ("correlation=", correlation)
    print ("friction=", friction)
    print ("params=", params, "(len(params)=%d" % len(params) )
    if params:
        if len(params) < 3 or len(params) > 4 :
            raise RuntimeError("CoolingsBCs: wrong number of params (%s)!!!" % str(params))
        
        pressure=float(params[0])
        dpressure=float(params[1])
        temperature=float(params[2])
        if len(params) == 4 :
            uguess=float(params[3])
            print ("CoolingBCs: pressure=%g, temperature=%g" % (pressure,temperature))
        else:
            uguess=0 
            print ("CoolingBCs: pressure=%g, temperature=%g" % (pressure,temperature))
    else:
        pressure=15 # bar
        dpressure=25 #bar
        temperature=30 # Celsius
        uguess=0
            
    print ("Pextra=", Pextra)
    print ("fguess=", fguess)
    print ("uguess=", uguess)

    # Set Therm Boundary conditions (!!!!Need to take care of units!!!!)
    CoolingBCs = {}
    (heat_coeff, Tw, flow) = self.cooling(Assembly, pressure, dpressure, temperature, correlation, friction, Pextra, fguess, uguess)

    print ("h: ", heat_coeff)
    print ("Tw: ", Tw)
    return CoolingBCs


#
# To operate from command line

if __name__ == "__main__":
    import os
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="name of the numerical model to be stored")
    parser.add_argument("--geometry", help="specify the geometry attached to the model", type=str)
    parser.add_argument("--inputcurrent", help="specify the input current (kA)", type=float)
    parser.add_argument("--CoolantTemp", help="specify the water input temperature (Celsius)", type=float, default=20)
    parser.add_argument("--dT", help="specify the dT (Celsius)", type=float, default=20) 
    parser.add_argument("--BP", help="specify the BP pressure (Bar)", type=float, default=10) 
    parser.add_argument("--dP", help="specify the dP (Bar)", type=float, default=15) 
    parser.add_argument("--correlation", help="specify the correlation used for heat exchange coefficient calculation (Colburn)", type=str, default="colburn") 
    parser.add_argument("--friction", help="specify the friction coefficient", type=float, default=0.055) 
    parser.add_argument("--headlosses", help="specify the extra pressure head losses", type=float, default=1) 
    parser.add_argument("--flow", help="specify the total flow rate", type=float) 
    parser.add_argument("--with_meanU", help="use mean water velocity", action="store_true")
    parser.add_argument("--verbose", help="activate verbose mode", action="store_true")
    args = parser.parse_args()

    meshdata = ""
    geometry = ""
    verbose = False

    if not args.geometry:
        geometry = args.name
    else:
        geometry = args.geometry
        
    if not args.meshdata:
        meshdata =  args.name
    else:
        meshdata = args.meshdata

    if args.verbose:
        verbose = args.verbose
        
    import getpass
    UserName = getpass.getuser()

    Assembly = Insert(geometry)
    Assembly.load()
    
    if verbose:
        print (Assembly)
    

    # Display Heat exchange Coefficient by Channel
    # if inputcurrent is given:
    #    compute Vpump
    #    compute DeltaT, Tw and DeltaP
    meanTw =  args.CoolantTemp + (args.dT/2.)
    meanP =  args.BP + (args.dP/2.)
    Pextra=args.headlosses
    fguess=args.friction

    # initial guess for h: Montgomery
    TotalCoolingSurf = 0
    maxDh = 0
    for i in range(len(Assembly.Helices)+1):
        (dh, L, S, Zmin, Zmax) = mymodel.hydraulic_diameter(i, Assembly)
        
        TotalCoolingSurf += S
        maxDh = max(dh, maxDh)
        
    meanU=0
    if not args.with_meanU:
        meanU = args.flow/TotalCoolingSurf
        meanH = 1426*(1+1.5e-2*meanTw)*math.pow(meanU,0.8)/math.pow(maxDh*1.e-3,0.2)
        print ("\nCst (Montgomery): meanU=%g, meanH=%g\n" % (meanU, meanH))
        
    
    for correl in ["montgomery", "colburn", "dittus", "silberberg"]:
        for friction in ["constant", "blasius", "filonenko","karman","rough"]:
            print ("#Channel\tDh[mm]\tL[mm]\tS[m²]\t\th[W/m2/K]\tu[m/s]\t\tRe[]\tPr[]\tCf[%s]\t%s" % (friction, correl) )
            print ("-----------------------------------------------------------------------------------")
            (h, T,  flow) = mymodel.cooling(Assembly, meanP, args.dP, meanTw, correl, friction, Pextra, fguess, uguess=meanU)
            print ("-----------------------------------------------------------------------------------")

            # loop over channels:
            Q = 0
            for i in range(len(h)):
                dh = (flow[i])[0]
                L = (flow[i])[1]
                S = (flow[i])[2]
                u = (flow[i])[3]
                Reynolds = (flow[i])[4]
                Prandtl = (flow[i])[5]
                Nusselt = (flow[i])[6]
                f = (flow[i])[7]
                print ("%d\t\t%g\t%g\t%g\t\t%g\t\t%g\t\t%g\t%g\t%g" % (i, dh, L, S, h[i], u, Reynolds, Prandtl, f))
                Q += S*u
            print ("FlowRate|l/s]: %g, Tw[C]: % g, P[Bar]: %g, extraP[Bar]: %g" % (Q/1.e+3, meanTw, meanP, Pextra) )
            print ("\n")
            
    
