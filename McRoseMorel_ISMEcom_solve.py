#---------------------------------------------------------------
# some needed python modules
#---------------------------------------------------------------
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------
# formation and dissociation rates for different model ligands
#---------------------------------------------------------------
kf_EDTA  = 7.2e4  # [1/(M hr)]
kf_EntB  = 3.6e9
kf_FeChr = 1.7e9
kf_FeOxA = 7.1e9
kd_EDTA  = 3.6e-3 # [1/hr]
kd_EntB  = 5.7e-2
kd_FeChr = 1.8e-4              # value from Witter et al., 2000
kd_FeChr_Boiteau = 3.6e-4      # value from Boiteau et al. 2022
kd_FeOxA = 5.4e-3              # value from Witter et al., 2000
kd_FeOxA_Boiteau = 1.3e-4      # value from Boiteau et al. 2022
FeEDTA0  = 100.0e-9
EntB0 = 50.0e-9

#---------------------------------------------------------------
# definition of the model equations
#---------------------------------------------------------------

# first model, assuming FeEDTA concentration is unchanged, so we solve
# only for the concentrations of Fe' and Fe bound to the added
# ligand. 'free' ligand is calculated from mass balance
def modelEntB(t, Y):
    Feprime    = Y[0]
    Fe_ligand  = Y[1]
    Lig_free = Lig_added - Fe_ligand
    if Lig_free < 0:
        Lig_free = 0
        Fe_ligand = Lig_added
    dFeprimedt = kd_EDTA * FeEDTA0 - kf_lig * Feprime * Lig_free + kd_lig * Fe_ligand
    dFeliganddt  = kf_lig * Feprime * Lig_free - kd_lig * Fe_ligand
    dYdt = [dFeprimedt, dFeliganddt]
    return dYdt

# second model, taking into account that FeEDTA decreases over time by dissociation
def model2EntB(t, Y):
    Feprime    = Y[0]
    Fe_ligand  = Y[1]
    FeEDTA     = Y[2]
    Lig_free = Lig_added - Fe_ligand
    if Lig_free < 0:
        Lig_free = 0.0
        Fe_ligand = Lig_added
    dFeprimedt = kd_EDTA * FeEDTA - kf_lig * Feprime * Lig_free + kd_lig * Fe_ligand
    dFeliganddt  = kf_lig * Feprime * Lig_free - kd_lig * Fe_ligand
    dFeEDTAdt  = -kd_EDTA * FeEDTA
    dYdt = [dFeprimedt, dFeliganddt, dFeEDTAdt]
    return dYdt

#---------------------------------------------------------------
# numerical solution of the model equations over 240 hours (10 days)
#---------------------------------------------------------------
if __name__ == "__main__":

    question = """
Which ligand is added? (Enter an integer number between 1 and 5)
   1: Enterobactin
   2: Ferrichrome (rates by Witter)
   3: Desferrioxamine B (rates by Witter)
   4: Ferrichrome (rates by Boiteau)
   5: Desferrioxamine B (rates by Boiteau)
"""
    answer = int(input(question))

    if (answer == 1):
        Lig_type = "Enterobactin"
        kf_lig = kf_EntB
        kd_lig = kd_EntB
    elif (answer == 2):
        Lig_type = "Ferrichrome (Witter)"
        kf_lig = kf_FeChr
        kd_lig = kd_FeChr
    elif (answer == 3):
        Lig_type = "Desferrioxamine B (Witter)"
        kf_lig = kf_FeOxA
        kd_lig = kd_FeOxA
    elif (answer == 4):
        Lig_type = "Ferrichrome (Boiteau)"
        kf_lig = kf_FeChr
        kd_lig = kd_FeChr_Boiteau
    elif (answer == 5):
        Lig_type = "Desferrioxamine B (Boiteau)"
        kf_lig = kf_FeOxA
        kd_lig = kd_FeOxA_Boiteau
    else:
        print("did not recognise input")
    print("Using formation and dissociation constants for " + Lig_type)
    
    # how much ligand are we adding?
    
    question = """
How much of the competing ligand is added?
(enter the concentration in nmol/L)
"""
    
    answer = float(input(question))
    Lig_added = answer * 1.0e-9 
    print("added ligand amount: ", Lig_added, " mol/L")
    
    tspan = [0, 240] # we integrate 10 days = 240 hours
    # over the first 2 hours we want output every 0.02 hours, after that every hour
    teval = np.concatenate( (np.arange(tspan[0], 2, 0.02), np.arange(2.0, tspan[1]+1, 1.0)) )
    
    # initial condition: Fe' and Fe bound to new ligand = 0
    InitCon = [0.0, 0.0]
    sol = solve_ivp(modelEntB, tspan, InitCon, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)
    
    # initial condition: Fe' and Fe bound to ligand added = 0, FeEDTA = 100 nM
    InitCon2 = [0.0, 0.0, FeEDTA0]
    sol2 = solve_ivp(model2EntB, tspan, InitCon2, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)
    
    #---------------------------------------------------------------
    # make a plot
    #---------------------------------------------------------------
    fig,ax = plt.subplots()
    ax.semilogy(sol.t, sol.y[0], label="Fe' (2 eqns)")
    ax.semilogy(sol.t, sol.y[1], label='FeEntB')
    ax.semilogy(sol2.t, sol2.y[0], '--', label="Fe' (3 eqns)")
    ax.semilogy(sol2.t, sol2.y[1], '--', label='FeEntB')
    ax.legend()
    ax.set_xlabel("t [hr]")
    ax.set_ylabel("Fe species [M]")
    plt.show()

