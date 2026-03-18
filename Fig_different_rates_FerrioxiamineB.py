#---------------------------------------------------------------
# some needed python modules
#---------------------------------------------------------------
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#---------------------------------------------------------------
# this is where the basic rate constants and model equations are defined
#---------------------------------------------------------------
import McRoseMorel_ISMEcom_solve as McR  

# how much ligand are we adding?

question = """
How much of the competing ligand is added?
(enter the concentration in nmol/L)
"""
    
answer = float(input(question))
McR.Lig_added = answer * 1.0e-9 
print("added ligand amount: ", McR.Lig_added, " mol/L")
    
tspan = [0, 240] # we integrate 10 days = 240 hours
# over the first 2 hours we want output every 0.02 hours, after that every hour
teval = np.concatenate( (np.arange(tspan[0], 2, 0.02), np.arange(2.0, tspan[1]+1, 1.0)) )

# now do two runs: one for constants for Ferrichrome by Witter, one by Boiteau
Lig_type = "Ferrioxiamine A (Witter)"
McR.kf_lig = McR.kf_FeOxA
McR.kd_lig = McR.kd_FeOxA_Witter
    
# initial condition: Fe' and Fe bound to ligand added = 0, FeEDTA = 100 nM
InitCon2 = [0.0, 0.0, McR.FeEDTA0]
sol1 = solve_ivp(McR.model2EntB, tspan, InitCon2, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)
    
Lig_type = "Ferrioxiamine (Boiteau)"
McR.kf_lig = McR.kf_FeOxA
McR.kd_lig = McR.kd_FeOxA_Boiteau
    
# initial condition: Fe' and Fe bound to ligand added = 0, FeEDTA = 100 nM
InitCon2 = [0.0, 0.0, McR.FeEDTA0]
sol2 = solve_ivp(McR.model2EntB, tspan, InitCon2, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)
    
#---------------------------------------------------------------
# make a plot
#---------------------------------------------------------------
fig,ax = plt.subplots()
ax.semilogy(sol1.t, sol1.y[0], label="Fe' (Witter)")
ax.semilogy(sol1.t, sol1.y[1], label='FeOxA')
ax.semilogy(sol2.t, sol2.y[0], '--', label="Fe' (Boiteau)")
ax.semilogy(sol2.t, sol2.y[1], '--', label='FeOxA')
ax.legend()
ax.set_xlabel("t [hr]")
ax.set_ylabel("Fe species [M]")
ax.grid()
plt.show()

