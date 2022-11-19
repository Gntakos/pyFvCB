import math

class FvCB:
    """
    The FvCB model of Photosynthesis, as being described in Bonen 2019:
    
    Leaf Photosynthesis (2019). 
    Climate Change and Terrestrial Ecosystem Modeling, 167-188
    
    https://doi.org/10.1017/9781107339217.012

    Author of this python version: Georgios Ntakos

    More info:
    ----------
    - The model includes product-limited assimilation (triose phosphate)
    - Calculation of electron transport rate (J) is based on the 
        smaller root of the quadratic equation equation [11.24]
    - Temprature correction is based on the Arrhenius function (with energy activation ΔHa)
    - There is an extra temperature correction for Vcmax, Jmax and Rd, based on the idea that
        thermal breakdown of biochemical processes occurs at higher temperatures (fH)
    """

    def __init__(self, Kc_25 = 404.9, Ko_25 = 278.4, extTp = None, sim_return = 1):
        """
        Initialization of the FvCB model.
        Initializes Kc, Ko and Tp with optional user input values or defaults.

        Assigns values to the activation energies, with suffix "_DHa" [J mol-1].
        Assigns ΔHd and ΔS for photosynthetic temperature inhibition.
        Assigns constants Oi and gas_constant

        Inputs [Optional]
        ------
        Kc_25:  Michaelis constant for CO2 [μmol mol-1]
        Ko_25:  Michaelis constant for 02 [mmol mol-1]
        extTp:  Rate of triose phosphate utilization [μmol m-2 s-1]
        sim_return: Output of the FvCB model ||1: only A, other number: all 5 parameters||

        Outputs
        -------
        []
        """

        self.Kc_DHa = 79430
        self.Ko_DHa = 36380
        self.Gx_DHa = 37830
        self.Vcmax_DHa = 65330
        self.Rd_DHa = 46390
        self.Jmax_DHa = 43540
        self.DHd = 150000
        self.DS = 490
        
        #Rubisco biochemistry parameters
        self.Kc_25 = Kc_25
        self.Ko_25 = Ko_25
        self.extTp = extTp
        
        #Constants
        self.gas_const = 8.314  # J K–1 mol–1
        self.Oi = 209 #intercellular concentration of O2 [mmol mol-1] 
                      #OR partial pressure of O2 [mbar] "not clear!"
        self.sim_return = sim_return #for output selection between only A or [A, Ac, Aj, Ap, Rd]
    
    def rubisco_limited(self, Vcmax, Ci, Gx, Kc, Ko):
        """
        The Rubisco-limited photosynthesis 
        (Dependence on maximum Rubisco activity through the Vcmax parameter)
        
        Inputs
        ------
        Vcmax: Maximum rate of carboxylation [μmol m-2 s-1]
        Ci: Intercellular CO2 [ppm]
        Gx: CO2 compensation point 'Γ*' [ppm]
        Kc: Michaelis constant for CO2 [μmol mol-1]
        Ko: Michaelis constant for 02 [mmol mol-1]

        Outputs
        ------
        Ac: Rubisco limited gross assimilation rate [μmol m-2 s-1]
        """

        num = Vcmax*(Ci - Gx)
        denom = Ci + Kc*(1 + self.Oi/ Ko)
        Ac = num/denom
        return Ac

    def electron_transport_rate(self, Jmax, par, theta):
        """
        The rate of electron transport is the smaller root of the quadratic equation:
        Θ * J^2 - (Iph2 + Jmax)* J + Ips2 * Jmax = 0

        Inputs
        ------
        Jmax: Maximum electron transport rate [μmol m-2 s-1]
        par: Photosynthetically active radiation [mmol m-2 s-1]
        theta: Curvature parameter, 'Θ' [-]

        Outputs
        ------
        J: Electron transport tate [μmol m-2 s-1]
        """

        ipar = 0.8*par # This is approximation of absorbed PAR
        phiPS2 = 0.7 # Quantum yield of photosystem II [mol mol–1]

        Ips2 = 0.5 * phiPS2 * ipar # Light utilized in electron transport by photosystem II [μmol m–2 s–1]
        num = (Ips2 + Jmax) - math.sqrt((Ips2 + Jmax)**2 - 4 * theta * Ips2 * Jmax)
        denom = 2*theta
        J = num/denom
        return J

    def light_limited(self, J, Ci, Gx):
        """
        Light limited photosynthesis. Or assimilation limited by RuBP regeneration

        Inputs
        ------
        J: Electron transport rate ||Through calling 'electron_transport_rate' method|| [μmol m-2 s-1]
        Ci: Intercellular CO2 [ppm]
        Gx: CO2 compensation point, 'Γ*' [ppm]

        Outputs
        ------
        Aj: Light limited gross assimilation rate [μmol m-2 s-1]
        """

        Aj = (J/4)*((Ci - Gx)/(Ci + 2*Gx))
        return Aj

    def arrhenius(self, temp, DHparam):
        """
        Arrhenius function for parameter temperature correction

        Inputs
        ------
        temp: Temperature [oC]
        DHparam: Activation energy of the corresponding parameter ||Assigned in __init__|| [J mol-1]

        Output
        ------
        t_corr: Temperature correction for corresponding parameter [-]
        """

        ktemp = temp + 273.15 #change to K for arrhenius temperature responce calculation
        t_corr = math.exp((DHparam / (298.15 * self.gas_const)) * (1 - (298.15 / ktemp)))
        return t_corr
    
    def hot_inhibition(self, temp):
        """
        Extra temperature correction for Vcmax, Jmax and Rd. This method is used to correct for 
        higher temepreatures

        Inputs
        ------
        temp: Temperature  [oC]

        Output
        ------
        hot_inhib: parameter for thermal breakdown of biochemical processes above 25 oC
        """

        ktemp = temp + 273.15 #change to K for arrhenius temperature responce calculation
        num = 1 + math.exp((298.15 * self.DS - self.DHd) / 298.15 * self.gas_const)
        denom = 1 + math.exp((ktemp * self.DS - self.DHd)/ ktemp * self.gas_const)
        hot_inhib = num/denom
        return hot_inhib
    
    def calcphot(self, CO2 = 266, PAR = 2000, temper = 25, Vcmax_25 = 60, Jmax_25 = 100.2, Gx_25 = 42.75, theta = 0.9, Rd_25 = 0.9):
        """
        "Main" method for calculation of Photosynthesis:
        1. Temperature correction of all the parameters (plus extra correction for Vcmax, Jmax & Rd)
        2. Calculation of rubisco limited assimilation rate
        3. Calculation of light limited assimilation rate
        4. Calculation of Triose phosphate based limited assimilation rate (product based)
        5. Calculate gross assimilation as the minimum of all the three rates above
        6. return 1 = A (Gross assimilation - Rd) or all the rates + Rd (A, Ac, Aj, Ap, Rd)

        Inputs
        ------
        CO2: CO2 concentration [ppm]
        PAR: Photosynthetically-active radiation [μmol m-2 s-1]
        temper: Temperature [oC]
        Vcmax_25: Maximum rate of carboxylation at 25 oC [μmol m-2 s-1]
        Jmax_25: Maximum electron transport rate at 25 oC [μmol m-2 s-1]
        Gx_25: CO2 compensation point, 'Γ*'at 25 oC [ppm]
        theta: Curvature parameter, 'Θ' [-]
        Rd_25: Dark respiration at 25 oC [μmol m-2 s-1]

        Output
        ------
        A: minimum of light, rubisco and product limited photosynthesis, minus dark respiration

        OR (based on the 'sim_return' parameter value)

        A, Ac, Aj, Ap, Rd
        """

        # Calculate temperature correction for all parameters
        Vcmax = Vcmax_25 * self.arrhenius(temper, self.Vcmax_DHa) * self.hot_inhibition(temper)
        Jmax = Jmax_25 * self.arrhenius(temper, self.Jmax_DHa) * self.hot_inhibition(temper)
        Rd = Rd_25 * self.arrhenius(temper, self.Rd_DHa) * self.hot_inhibition(temper)
        Kc = self.Kc_25 * self.arrhenius(temper, self.Kc_DHa)
        Ko = self.Ko_25 * self.arrhenius(temper, self.Ko_DHa)
        Gx = Gx_25 * self.arrhenius(temper, self.Gx_DHa)
        
        # Calculate Assimilation Responses for ligh, CO2, and triose phosphate
        Ac = self.rubisco_limited(Vcmax, CO2, Gx, Kc, Ko)
        J = self.electron_transport_rate(Jmax, PAR, theta)
        Aj = self.light_limited(J, CO2, Gx)
        
        if self.extTp is None:
            Tp = 0.167 * Vcmax
        else:
            Tp = self.extTp

        Ap = 3 * Tp
        
        GASS = min(min(Ac, Aj), Ap)
        A = GASS - Rd
        #return Ac, Aj, Ap, A
        return A if self.sim_return == 1 else A, Ac, Aj, Ap, Rd




