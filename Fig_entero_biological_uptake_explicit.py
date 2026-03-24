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

# which ligand we are using
Lig_type = "Enterobactin"
McR.kf_lig = McR.kf_EntB
McR.kd_lig = McR.kd_EntB

# how sensitive is SAR 11 growth to iron (kFe for Fe')
McR.kfe_bac = 5e-12
McR.mu_bac  = 0.54/24

# now do three model integrations, with three different amounts of
# added ligand

tspan = [0, 240] # we integrate 10 days = 240 hours
# over the first 2 hours we want output every 0.02 hours, after that every hour
teval = np.concatenate( (np.arange(tspan[0], 2, 0.02), np.arange(2.0, tspan[1]+1, 1.0)) )

# initial cell number
initial_cell = 2.0e-7

# first with 25 nM added Enterobactin
McR.Lig_added = 25.0e-9 
InitCon2 = [0.0, 0.0, McR.FeEDTA0, McR.Lig_added, initial_cell]
sol1 = solve_ivp(McR.model4EntB, tspan, InitCon2, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)

# then with 50 nM added Enterobactin
McR.Lig_added = 50.0e-9 
InitCon2 = [0.0, 0.0, McR.FeEDTA0, McR.Lig_added, initial_cell]
sol2 = solve_ivp(McR.model4EntB, tspan, InitCon2, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)

# finally with 100 nM added Enterobactin
McR.Lig_added = 100.0e-9 
InitCon2 = [0.0, 0.0, McR.FeEDTA0, McR.Lig_added, initial_cell]
sol3 = solve_ivp(McR.model4EntB, tspan, InitCon2, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)

#---------------------------------------------------------------
# calculate cumulative biological Fe uptake
#---------------------------------------------------------------
#C0 = 2.0e7      # initial cell number of SAR-11 in cells/L
#mu = 0.5/24.0   # unlimited growth rate of SAR-11 in the experiments in 1/h
#ucell = 2.4e-21 # cellular uptake rate in mol Fe/cell/h
#BFe = C0 * ucell / mu * np.exp( teval * mu )

#---------------------------------------------------------------
# make a plot
#---------------------------------------------------------------
fig,ax = plt.subplots(2,1)

ax[0].semilogy(sol1.t, sol1.y[0], 'b-', label="Fe' with 25 nM EntB")
ax[0].semilogy(sol2.t, sol2.y[0], 'b--', label="Fe' with 50 nM EntB")
ax[0].semilogy(sol3.t, sol3.y[0], 'b-.', label="Fe' with 100 nM EntB")
# ax.semilogy(teval, BFe, 'r-', label="cumulative Fe uptake")
ax[0].legend()
ax[0].set_xlabel("t [hr]")
ax[0].set_ylabel("Fe species [M]")
ax[0].grid()

ax[1].semilogy(sol1.t, sol1.y[4], 'b-', label="with 25 nM EntB")
ax[1].semilogy(sol2.t, sol2.y[4], 'b--', label="with 50 nM EntB")
ax[1].semilogy(sol3.t, sol3.y[4], 'b-.', label="with 100 nM EntB")
# ax.semilogy(teval, BFe, 'r-', label="cumulative Fe uptake")
ax[1].legend()
ax[1].set_xlabel("t [hr]")
ax[1].set_ylabel("cell numbers [#/L]")
ax[1].grid()

plt.savefig("uptake_explicit_enterobactin_5pM.png")
plt.show()
