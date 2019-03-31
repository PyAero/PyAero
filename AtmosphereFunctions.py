# -*- coding: utf-8 -*-
"""
Standard Athmosphere


"""

def StandardAtmosphere(state,Altitude):

    #Constants    
    h       = Altitude/3.2808399
    C2K     = 273.15
    e       = 2.71828

    #Troposphere
    if h < 11000:
        
        T0      = 15.04     #C intercept Temp at 0m
        dtdh    = -0.00649  #dC/dm 
        
        T = T0 + dtdh*h                 #C
        
        Tref    = 288.08    #K ref Temp
        Pref    = 101.29    #KPa ref Pressure
        eT      = 5.256     #  Temperature ratio pressure exponent
        
        P = Pref*((T+C2K)/Tref)**eT     #KPa
    
    #Lower Stratosphere
    elif h < 25000:
        
        T = -56.46                      #C Temp for 11000 to 25000
        
        h0      = 1.73      #m ref Altitude
        dtdh    = 0.000157  #dC/dm
        Pref    = 22.56     #KPa ref Pressure
        
        P = Pref*e**(h0-dtdh*h)         #KPa
    
    #Upper Stratosphere    
    elif h < 40000:
    
        T0      = -131.21   #C intercept Temp at 0m
        dtdh    = 0.00299   #dC/dm 
        
        T = T0 + dtdh*h                 #C
        
        Tref    = 216.6     #K ref Temp
        Pref    = 2.488     #KPa ref Pressure
        eT      = -11.388   #  Temperature ratio pressure exponent
        
        P = Pref*((T+C2K)/Tref)**eT     #KPa
    
    #Catch
    else:
        
        return "Function only goes to 120 kft"
        
    # Density 
    Rair        = 0.2869    #KPa-m^3/kg-K - Air Gas Constant
    
    rho = P/(Rair*(T + C2K))          #kg/m^3  

    #Speed of Sound
    gamma = 1.4
    a = (gamma*1000*P/rho)**0.5  #m/s
    
    
    P   = P/6.89476         #psi
    T   = 9/5*T + 32        #F
    rho = rho/515.379       #slug/ft^3
    a   = a*3.28084         #ft/s
    
    #Output
    if   state == "Pres":
        return P
    elif state == "Temp":
        return T
    elif state == "Rho":
        return rho
    elif state == "a":
        return a
    else:
        return [P,T,rho,a]

    
def MachAltKeas(Alt,Keas):
    
   rho0     = StandardAtmosphere("Rho",0)
   rho      = StandardAtmosphere("Rho",Alt)    
   Va       = StandardAtmosphere("a",Alt)
    
   Ktas     = Keas*(rho/rho0)**0.5
    
   Mach     = Ktas/Va
    
   return Mach


def KtasMachAlt(Mach,Alt):

   Va       = StandardAtmosphere("a",Alt)
   
   Ktas     = Mach*Va
   
   return Ktas
   
def KeasMachAlt(Mach,Alt):

   Va       = StandardAtmosphere("a",Alt)
   Ktas     = Mach*Va
   
   rho0     = StandardAtmosphere("Rho",0)
   rho      = StandardAtmosphere("Rho",Alt)   
   
   Keas     = Ktas*(rho/rho0)**0.5
   
   return Keas

def AltKeasMach(Keas,Mach):
        
    a1 = 0 #ref alt 1
    a2 = 10000 #ref alt 2
    
    k1 = KeasMachAlt(Mach,a1)
    k2 = KeasMachAlt(Mach,a2)
    
    ai = (a2-a1)/(k2-k1)*(Keas-k1) + a1 #Newton-Rhapson
         
    ki = KeasMachAlt(Mach,ai)
    
    error = (Keas-ki)/Keas
            
    while abs(error) > 0.00001:   #Drive down error
    
        a2 = a1
        a1 = ai
    
        k1 = KeasMachAlt(Mach,a1)
        k2 = KeasMachAlt(Mach,a2)
        
        ai = (a2-a1)/(k2-k1)*(Keas-k1) + a1 #Newton-Rhapson
             
        ki = KeasMachAlt(Mach,ai)
        
        error = (Keas-ki)/Keas    
    
    
    Alt = ai

    return Alt    
    
def qMachAlt(Mach,Alt):

    Ktas = KtasMachAlt(Mach,Alt)
    rho  = StandardAtmosphere("Rho",Alt)
    
    q = 0.5*rho*Ktas**2
    
    q = q/144
    
    return q

def P0MachAlt(Mach,Alt):
    
    p   = StandardAtmosphere("Pres",Alt)
    gamma = 1.4
    
    a   = (gamma-1)/2
    b   = gamma/(gamma-1)
    
    P0  = p*(1+a*Mach**2)**b
    
    return P0