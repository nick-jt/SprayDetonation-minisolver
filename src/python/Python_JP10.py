#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cantera as ct
import numpy as np
import math
import sys
from scipy.integrate import *
import matplotlib.pyplot as plt
from sdtoolbox.thermo import soundspeed_fr
import signal


# # 1-D Steady Spray Detonation Solver

# ## This data structure is used to solve our system:
# * **__init__()** is an operator which initializes the constant parameters of our system.
# * **getLatentHeat(Td,w_f)** is a function which gets the latent heat of droplet vaporization
# * **getVaporPressure(Td)** get's the vapor pressure using the user's preferred method.
# * **StateVectorFunction(x,y)** is the system of ODEs to be integrated. State terms are inputted, and this function evaluates the evolution terms. The python solve_ivp function stored in the integrator function will then integrate this.
# * **integrator(D)** will generate the initial conditions using the parameters provided, solve the postshock conditions, then set up and solve our IVP.
# 
# Additional test section modifications have been made. All modified vaporization functions are added through inheritance & method override.

# In[2]:


class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class VaporizationError(Error):
    """Error associated with non-physical vaporization calculations."""
    def __init__(self,message):
        self.message = message

class ParameterError(Error):
    def __init(self,message):
        self.message = message

class TimeoutError(Exception):
    """Raised when our runtime exceeds a value"""
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutError("Function timed out...")

class Detonation:

    def __init__(self):
        self.lchar  = 3.81 * 0.01 / 4 # Characteristic Length (m)
        self.Pr     = 1                 # Prandtl Number # TODO: Consider changing to 0.67
        self.Le     = 1                   # Lewis Number
        self.Tw     = 298              # Wall Temperature (K)
        self.rhod   = 285             # Density of droplets Engineering toolbox
        self.T_in   = 298            # Preshock temperature (K)
        self.P_in   = 1e5          # Preshock pressure (Pa)
        self.Length = 1.0            # Domain length
        self.oxidizer = "O2:0.21145,N2:0.78855"# air composition
        
    def getLatentHeat(self,Td,w_f,Pfs): 
        """
            can we just calculate from vapor pressure?
            "The enthalpies of vaporization and sublimation
            of exo- and endo-tetrahydrodicyclopentadienes at 
            T=298.15K" - University of Missouri and NASA Glenn
        """
        return 541911.764706
        
    def getVaporPressure(self,Td):
        """
            "HIGH-TEMPERATURE HEAT-CAPACITY MEASUREMENTS AND CRITICAL
            PROPERTY DETERMINATIONS USING A DIFFERENTIAL SCANNING CALORIMETER"
            WPAFB
        """
        A = 2.39742
        B = -1.28961
        C = 0.90779
        Tc = 698
        Pc = 3733000
        Tr = Td/Tc
        P = Pc*np.exp((1-1/Tr)*np.exp(A+B*Tr+C*Tr**2))
        return P
    
    def getSpecificHeat(self, T, Rsp):
        c_0 = 3.3218
        c_1 = 0.07975
        c_2 = 27.6975
        c_3 = 1470
        Cpd = Rsp*(c_0 + c_1*T**0.85 + c_2*(c_3/T)**2              * np.exp(c_3/T)/(np.exp(c_3/T)-1)**2 )
        return Cpd
    
    def postshock(self, V, P_in, T_in, q, mech):
        gas = ct.Solution(mech)
        gas.TPX = T_in,P_in,q
        p1 = gas.P
        rho1 = gas.density
        u1 = V
        h1 = gas.enthalpy_mass
        
        rho2 = 5 # initial guess
        rho2_prev = 4
        
        while np.abs(rho2-rho2_prev)>1e-4*rho2:
            u2 = rho1*u1/rho2
            p2 = p1+rho1*u1**2-rho2*u2**2
            h2 = h1+0.5*u1**2-0.5*u2**2
            gas.HP = h2,p2
            rho2_prev = rho2
            rho2 = gas.density
            
        return gas
    
    def getDragCoefficient(self,Ma,Re): 
        """Compressibility and Rarefaction Effects on Drag of a Spherical Particle"""
        if Re < 0.1:
            Cdd = 24/Re
        elif Re < 45:
            Cdd = (24/Re)*(1+0.15*Re**0.687)
        elif Re > 45:
            if Ma < 0.89:
                Gm = 1-1.525*Ma**4
            elif Ma >= 0.89:
                Gm = 10**-4 * (2 + 8*np.tanh( 12.77*(Ma-2.02) ))
            if Ma < 1.45:
                Cm = 1/3 * (5 + 2*np.tanh(3*np.log(Ma+0.1)))
            elif Ma > 1.45:
                Cm = 2.044 + 0.2*np.exp(-1.8*np.log(Ma/1.5)**2)
            Hm = 1-0.258*Cm/(1+514*Gm)
            Cdd = Hm * (24/Re)*(1+0.15*Re**0.687) + 0.42*Cm/(1+42500*Gm*Re**(-1.16))
        return Cdd

    def ErrMngmt(self,y):
        for P,i in enumerate(y):
            if P<0:
                raise ParameterError("ERROR: y["+str(i)+"] is negative at x="+str(x))
        return

        
    def StateVectorFunction(self,x,y):

        """
        This function is called by the integrator to determine evolution terms.
        
        INPUTS:
            x = distance
            y = [time, gas velocity,droplet velocity, gas temperature, gas density,
                droplet radius, droplet temperature, gas species mass fractions]
        OUTPUTS:
            Numpy array of state evolution terms.
        """

        self.ErrMngmt(y)

    #    UNPACKING VARIABLES
        time        = y[0]
        ug          = y[1]
        ud          = y[2]
        Tg          = y[3]
        rhog        = y[4]
        rd          = y[5]
        Td          = y[6]
        Yg          = y[7:]
        
    #    ACQUIRING DATA FROM GLOBAL
        Cdw         = self.Cdw
        Chw         = self.Chw
        Rd0         = self.Rd
        lchar       = self.lchar
        Pr          = self.Pr
        Le          = self.Le
        Tw          = self.Tw
        rhod        = self.rhod
        nu0         = self.nu0
        D           = self.V
        lam         = self.lam
        alpha       = self.alpha
        gas         = self.gas
        gas_Td      = self.gas_Td
        fuel_index  = self.fuel_index

        
    #    FUNDAMENTAL PARAMETERS
        gas.TDY     = Tg,rhog,Yg
        nd          = nu0 / ud
        gam         = gas.cp_mass / gas.cv_mass
        grs         = gam / (gam - 1)
        w           = gas.mean_molecular_weight
        Rsp         = ct.gas_constant/w
        w_k         = gas.molecular_weights
        w_f         = w_k[fuel_index]
        omega       = gas.net_production_rates
        M           = ug / np.sqrt(gam*ct.gas_constant/w*gas.T)
        droplet_y   = np.linspace(0,0,gas.n_species)
        droplet_y[fuel_index] = 1
        
        if rd > 1e-2 * Rd0 and alpha > 0 and nd > 1e-10*100**3:
        
        #    DROPLET RELATED PARAMETERS
            Pfs = self.getVaporPressure(Td)
            L = self.getLatentHeat(Td,w_f,Pfs)
            Cpd = self.getSpecificHeat(Td,Rsp)
            Xfs = Pfs / gas.P
            if L>1e6 or L<1e5:
                raise VaporizationError("Latent heat fell outside range: (100,000 to 1,000,000 J/kg). L="+str(L))
            if Xfs>1 or Xfs<0:
                raise VaporizationError("Surface fuel mole fraction falls outside of range: "+str(Xfs))
                
            W_nofuel = (1 - gas.Y[fuel_index]) / (1 / w - gas.Y[fuel_index] / w_f)
            Yfs = Xfs * w_f / (Xfs * w_f + (1 - Xfs) * W_nofuel)
            By = (Yfs - gas.Y[fuel_index]) / (1 - Yfs)
            Bh = gas.cp_mass * (Tg - Td) / L

        #    DROPLET HEATING AND VAPORIZATION WITH CORRECTION
            Re = rhog * np.abs(ud - ug) * 2 * rd / gas.viscosity
            mdotv = nd * 4 * np.pi * rd * lam / (Le * gas.cp_mass) *                  np.log(1 + By) * (1 + 0.276 * Re**0.5 * Pr**(1/3))
            qd = nd * 4 * np.pi * rd * lam / gas.cp_mass * np.log(1 + Bh) *                 L * (1 + 0.276 * Re**0.5 * Pr**(1/3))
            
            # Random tests
            mdotv = mdotv*1.0
            qd = qd*1.0
            
        #    LOSSES
            Cdd = 22 * Re**(-1) * (1 + 0.276 * Re**0.5 * Pr**(1/3))
            fd = nd * Cdd * 4 * np.pi * rd**2 * rhog * np.abs(ud - ug) * (ud - ug) / 2
            #Mach_dr = np.abs(ud-ug) / np.sqrt(gam*ct.gas_constant/w*gas.T) # Droplet mach number
            #Cdd = self.getDragCoefficient(Mach_dr,Re)
            #fd = nd * Cdd * np.pi * rd**2 * rhog * np.abs(ud - ug) * (ud - ug) / 2 # Removed a factor of four
            
        #    DROPLET EVOLUTION
            drddx = -mdotv / (rhod * 4 * np.pi * rd**2 * nu0)
            duddx = -fd / (rhod * nd * 4/3 * np.pi * rd**3 * ud)
            dTddx = (qd - mdotv * L) / (rhod * nd * 4/3 * np.pi * rd**3 * ud * Cpd)
            gas_Td.TDY = Td,rhog,Yg
            enth = (gas.partial_molar_enthalpies[fuel_index] -                 gas_Td.partial_molar_enthalpies[fuel_index]) / w_f
            
        else:
                
            rd = ud = mdotv = qd = drddx = duddx = dTddx = fd = enth = 0

    #    EXTERNAL LOSSES
        fw = Cdw / lchar * rhog * np.abs(D - ug) * (D - ug) / 2
        qw = Chw / lchar * rhog * np.abs(D - ug) * gas.cp_mass * (Tg - Tw)

    #    GAS VELOCITY EVOLUTION
        hgk = gas.partial_molar_enthalpies / w_k
        reaction_source = sum( (hgk - gas.cp_mass*Tg*w/w_k) *omega*w_k)
        S = fw * (grs*ug-D) + qw + fd*(grs*ug-ud) + qd + reaction_source +              mdotv * (grs*ug*(ud-ug) - gas.cp_mass*Tg*w/w_f -              (ud**2 - ug**2) / 2 + enth)
        dugdx = (gam - 1) * M**2 * S / ((M**2 - 1) * rhog * ug**2)

    #    SPECIES EVOLUTION
        dygdx = 1 / (rhog * ug) * (omega * w_k + mdotv * (droplet_y - gas.Y))
                    
    #    GAS TEMPERATURE EVOLUTION
        dTgdx = -ug / gas.cp_mass * dugdx + 1 / (rhog * ug * gas.cp_mass) *             (fd * ud + fw * D - sum(hgk * omega * w_k) - qw - qd             + mdotv * ((ud**2 - ug**2) / 2 - enth))

    #    GAS DENSITY EVOLUTION
        drhogdx = -rhog / ug * dugdx + mdotv / ug
        
        if ud >0 :
            dtddx = 1/ud
        else:
            dtddx = 0
        
        out = np.hstack((1/ug,dugdx,duddx,dTgdx,drhogdx,drddx,dTddx,dygdx))
        return out
        
    def Choked(self,x,y):
        self.gas.TDY = y[3],y[4],y[7:]
        gam = self.gas.cp_mass/self.gas.cv_mass
        w = self.gas.mean_molecular_weight
        return y[1]/np.sqrt(gam*ct.gas_constant/w*y[3])-0.999
    Choked.terminal = True
    Choked.direction = 1
        
    def integrator(self,V,myrtol=1e-6,myatol=1e-6):
        """
            This function first gets the postshock conditions, then integrates the IVP.
            
            Order: time, Ug, Ud, Tg, rhog, Rd, Td, Y
        """
        
        x=self.gas.n_atoms(self.fuel, 'C')
        y=self.gas.n_atoms(self.fuel, 'H')
        a = x + y/4

        # Initializing actual gas
        self.fuel_index   = self.gas.species_index(self.fuel)
        self.gas.TPX = self.T_in,self.P_in,self.fuel+':'+str(self.phi*(1-self.alpha))+', O2:'+str(a)+', N2:'+str(a*3.76)
        self.q       = self.gas.X
        self.rho0    = self.gas.density
        
        # Initializing fully vaporized gas
        self.gas_withfuel.TPX = self.T_in,self.P_in,self.fuel+':'+str(self.phi)+', O2:'+str(a)+', N2:'+str(a*3.76)
        self.q_withfuel = self.gas_withfuel.X
        
        # INITIAL DROPLET NUMBER DENSITY
        if self.Rd:
            MF = self.gas_withfuel.Y[self.fuel_index]*self.alpha / (self.gas_withfuel.Y[self.fuel_index]*(1-self.alpha) + sum(self.gas_withfuel.Y[0:self.fuel_index]) + sum(self.gas_withfuel.Y[self.fuel_index+1:]))
            self.nd0 = self.gas.density*MF/(self.rhod*4/3*np.pi*self.Rd**3)
        else: self.nd0 = 0
        self.nu0    = self.nd0*V
        self.V      = V
        
        # DETERMINE POSTSHOCK STATE
        self.gas = self.postshock(V, self.P_in, self.T_in, self.q, self.mech)
        self.lam = self.gas.thermal_conductivity
        u = self.rho0 / self.gas.density * V
        
        # SETTING UP INTEGRATOR
        Range   = [0,self.Length]
        IC      = np.hstack((0,u,V,self.gas.T,self.gas.density,self.Rd,self.T_in,self.gas.Y))
        self.sol = solve_ivp(self.StateVectorFunction,Range,IC,events=(self.Choked),method='LSODA',rtol=myrtol,atol=myatol,max_step=self.dx,min_step=self.dxmin)
        return self.sol
        
    def getSSvelocity(self,lower,upper,myrtol=1e-6,myatol=1e-6,itertol=1e-1):
        iter = 0
        V = (lower+upper)/2
        Vprev = V + 100
        FullDomain = False
        Error = False
        signal.signal(signal.SIGALRM, timeout_handler) # Handling timing
        while (np.abs(Vprev-V)>itertol or FullDomain==False) and Error==False:
            iter += 1
            print("     Iteration=%d, V=%.2f m/s"%(iter,V),flush=True)
            perturbcount = 0
            while True:
                signal.alarm(30) # set a timer for 60 seconds
                try:
                    sol = self.integrator(V,myrtol=myrtol,myatol=myatol)
                except (ct.CanteraError, VaporizationError) as err: # Sometimes the integrator fails
                    perturbcount +=1 
                    if perturbcount > 5:
                        print("ERROR")
                        Error = True
                        break
                    print("Perturbing: {0}".format(err),flush=True)
                    V += 1e-2*(V-Vprev)
                    continue
                except (TimeoutError) as err:
                    perturbcount += 1
                    if perturbcount > 5:
                        print("ERROR")
                        Error = True
                        break
                    print("Perturbing due to TimeoutError: {0}".format(err),flush=True)
                    V += 1e-2*(V-Vprev)
                    continue
                else:
                        signal.alarm(0)
                break
            if self.sol.t[-1]>=self.Length:
                upper = V
                FullDomain = True
            else:
                lower = V
                FullDomain = False
            Vprev = V
            V = (lower+upper)/2
        return (lower+upper)/2,sol

    def bracketMethodCDW(self,V,CDWl,CDWu,myrtol=1e-6,myatol=1e-6): # get SS drag coefficient given velocity
        self.V = V
        CDWmean = (CDWl+CDWu)/2
        FullDomain = False
        Error = False
        signal.signal(signal.SIGALRM, timeout_handler) # Handling timing
        while (CDWu-CDWl>1e-3 or FullDomain==False) and Error==False:
                print("%f %f"%(CDWmean,V),flush=True)
                self.Cdw = CDWmean
                nperturbs = 0
                while True:
                        signal.alarm(600) # set a timer for 60 seconds
                        try:
                                sol = self.integrator(self.V,myrtol=myrtol,myatol=myatol)
                        except ct.CanteraError:
                                nperturbs += 1
                                if nperturbs > 5:
                                        print("ERROR number of perturbations exceeds 5")
                                        Error = True
                                        break
                                print("Perturbing due to CanteraError",flush=True)
                                self.Cdw += 1e-2*(CDWmean-CDWl)
                                continue
                        except TimeoutError:
                                nperturbs += 1
                                if nperturbs > 5:
                                        print("ERROR number of perturbations exceeds 5")
                                        Error = True
                                        break
                                print("Perturbing due to TimeoutError",flush=True)
                                self.Cdw += 1e-2*(CDWmean-CDWl)
                                continue
                        else:
                                signal.alarm(0)
                        break
                if self.sol.t[-1]>=self.Length:
                        CDWu = CDWmean
                        FullDomain = True
                else:
                        CDWl = CDWmean
                        FullDomain = False
                CDWmean = (CDWl+CDWu)/2
        return  CDWmean,sol
    
    def SSvel_vs_CDW(self,Vlist,CDWlow,CDWhigh,myrtol=1e-6,myatol=1e-6):
        CDWlist = []
        for i,V in enumerate(Vlist):
            print("\nV=%f"%V)
            CDWmean,sol = self.bracketMethodCDW(V,CDWlow,CDWhigh,myrtol,myatol)
            CDWlist.append(CDWmean)
            print("Cdw=%.10f"%CDWmean)
        return CDWlist
    
    def CDW_vs_SSvel(self,CDWlist,Vlow,Vhigh,myrtol=1e-6,myatol=1e-6):
        Vlist = []
        V = Vhigh
        for i,CDW in enumerate(CDWlist):
            print("\nCDW=%f"%CDW)
            D.Cdw = CDW
            V,sol = self.getSSvelocity(Vlow,V,myrtol=1e-6,myatol=1e-6,itertol=1e-1)
            Vlist.append(V)
            print("\nV=%f"%V)
        return Vlist
    
    def plot(self,axis,var,fs=8,lbl="",mysol=""):
        xlim= self.xlim
        xscale = self.xscale
        nsp = self.gas.n_species
        idx_options = { 'Ug':1,'Ud':2,'Tg':3,'rhog':4,'Rd':5,'Td':6,'Yf':7+self.gas.species_index(self.fuel),                            'YO':7+self.gas.species_index('O2')}
        label_options = {   'Ug'    : 'Gas Velocity (m/s)',
                            'Ud'    : 'Droplet Velocity (m/s)',
                            'Tg'    : 'Gas Temperature (K)',
                            'rhog'  : 'Gas Density (kg/m^3)',
                            'Rd'    : 'Droplet Radius (m)',
                            'Rd^2'  : 'R_d ^ 2 (m^2)',
                            'Td'    : 'Droplet Temperature (K)',
                            'Yf'    :  str(self.fuel)+' Mass Fraction',
                            'YO'    : 'Oxygen Mass Fraction',
                            'M'     : 'Mach Number',
                            'HRR'   : 'HRR (W/m^3)',
                            'PHI'   : 'Equivalence Ratio',
                            'Therm' : 'Thermicity',
                            'P'     : 'Pressure (Pa)',
                            'Dd'    : 'Droplet Diameter (m)',
                            'Dd^2'  : 'D_d ^ 2 (m^2)',
                            'Xd'    : 'Droplet Distance from Shock (m)'}
        
        idx = idx_options.get(var)
        label = label_options.get(var)
        axis.set_xlim(xlim[0],xlim[1])
        
        if idx != None:
            data = mysol.y[idx]
        else:
            if var=='M':
                Marray = []
                for i in range(len(mysol.t)):
                    self.gas.TDY=mysol.y[3,i],mysol.y[4,i],mysol.y[7:,i]
                    c = np.sqrt(self.gas.cp_mass/self.gas.cv_mass*                                ct.gas_constant/self.gas.mean_molecular_weight*self.gas.T)
                    Marray.append(mysol.y[1,i]/c)
                data = Marray
            elif var=='HRR':
                HRRarray = []
                for i in range(len(mysol.t)):
                    self.gas.TDY=mysol.y[3,i],mysol.y[4,i],mysol.y[7:,i]
                    #HRRarray.append(-np.dot(self.gas.net_production_rates,self.gas.partial_molar_enthalpies))
                    HRRarray.append(-np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy))
                data = HRRarray
                intHR = trapz(HRRarray,mysol.y[0,:])
                print("ODE Integrated HR = %.10f J/m^3"%intHR)
            elif var=='Therm':
                thermicityArray = []
                
                for i in range(len(mysol.t)):
                    self.gas.TDY=mysol.y[3,i],mysol.y[4,i],mysol.y[7:,i]
                    w = self.gas.molecular_weights
                    hs = self.gas.standard_enthalpies_RT*ct.gas_constant*self.gas.T/w
                    dydt = self.gas.net_production_rates*w/self.gas.density

                    Therm = sum((self.gas.mean_molecular_weight/w
                                      -hs/(self.gas.cp_mass*self.gas.T))*dydt)
                    thermicityArray.append(Therm)
                data = thermicityArray
            elif var=='PHI':
                PHIarray = []
                for i in range(len(mysol.t)):
                    self.gas.TDY=mysol.y[3,i],mysol.y[4,i],mysol.y[7:,i]
                    gas = self.gas
                    phi = sum([gas.X[i]*(2*gas.n_atoms(i,'C')+0.5*gas.n_atoms(i,'H')) for i in range(gas.n_species)])                             / sum([gas.X[i]*gas.n_atoms(i,'O') for i in range(gas.n_species)])
                    PHIarray.append(phi)
                data = PHIarray
            elif var=='P':
                Parray = []
                for i in range(len(mysol.t)):
                    self.gas.TDY=mysol.y[3,i],mysol.y[4,i],mysol.y[7:,i]
                    P = self.gas.P
                    Parray.append(P)
                data = Parray
            elif var=='Dd':
                data = mysol.y[5,:]*2
            elif var=='Dd^2':
                data = (mysol.y[5,:]*2)**2
            elif var=='Rd^2':
                data = (mysol.y[5,:])**2
            elif var=='Xd':
                data = cumtrapz(mysol.y[2,:],mysol.y[0,:])
                data = np.append(data,data[-1])+2.41597065e-6
            else:
                print("Variable not found")
        
        
        figname = var+'.png'
        
        axis.set_xscale(xscale)
        axis.plot(mysol.t,data,label=lbl)
        axis.set_ylabel(label,fontsize=fs)
        axis.set_xlabel('Distance from shock (m)',fontsize=fs)
#         axis.set_xlabel('Time after shock (s)',fontsize=fs)
        axis.tick_params(axis='both', which='major', labelsize=fs)


# ## Initializations

# In[3]:

if __name__ == "__main__":

   i = int(sys.argv[1])

   # Initializing Detonation instance
   D = Detonation()
   D.mech  = 'JP10skeletal.yaml' # 40 species skeletal
   D.fuel  = 'C10H16'
   
   # Declaring initial gas state
   D.gas     = ct.Solution(D.mech)
   D.gas_Td  = ct.Solution(D.mech)
   D.gas_withfuel = ct.Solution(D.mech)
   
   
   # ## Liquid Detonation
   
   # In[4]:
   
   
   D.Chw = 0
   D.Cdw = 0
   D.phi = 1.0
   D.alpha = 1
   D.Dd = 5e-6
   D.Rd = D.Dd/2
   D.dx = 8.1064e-05 # Max_step
   D.dxmin = 1e-13
   # sol = D.integrator(1669.5078360997054,myrtol=1e-10,myatol=1e-15)
   # V,sol = D.getSSvelocity(1600,1840)#,myrtol=1e-12,myatol=1e-15)
   # CDW,sol = D.bracketMethodCDW(1300,0,0.02)
   # Vlist = np.linspace(1200,1800,10)
   # CDWlist = D.SSvel_vs_CDW(Vlist,0,0.06,myrtol=1e-4,myatol=1e-6)
   
   vels = np.linspace(1200,1800,100)
   V = vels[i]
   print(V)
   CDW,sol = D.bracketMethodCDW(V,0,0.20)
   print(CDW,V)
   
