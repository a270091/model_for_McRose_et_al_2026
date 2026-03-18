from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import McRoseMorel_ISMEcom_solve as McR

# which ligand we are using? 
Lig_type = "Enterobactin"
McR.kf_lig = McR.kf_EntB
McR.kd_lig = McR.kd_EntB
McR.kdeg_lig = 1.0/(10*24)

tspan = [0, 240] # we integrate 10 days = 240 hours
# over the first 2 hours we want output every 0.02 hours, after that every hour
teval = np.concatenate( (np.arange(tspan[0], 2, 0.02), np.arange(2.0, tspan[1]+1, 1.0)) )

# run with 4 equations
McR.Lig_added = 50.0e-9 
InitCon = [0.0, 0.0, McR.FeEDTA0, McR.Lig_added]
sol1 = solve_ivp(McR.model3EntB, tspan, InitCon, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)

# run with 3 equations
McR.Lig_added = 50.0e-9 
InitCon2 = [0.0, 0.0, McR.FeEDTA0]
sol2 = solve_ivp(McR.model2EntB, tspan, InitCon2, method='BDF', rtol=1.0e-12, atol=1.0e-20, t_eval=teval)

fig,ax = plt.subplots()
ax.plot(sol1.t, sol1.y[0], label="Fe' with apo-Ent decay")
ax.plot(sol1.t, sol1.y[1], label='FeEnt')
ax.plot(sol1.t, sol1.y[3], label='apo-Ent')
ax.plot(sol2.t, sol2.y[0], '--', label="Fe' without decay")
ax.plot(sol2.t, sol2.y[1], '--', label='FeEnt')
ax.legend()
ax.set_xlabel("t [hr]")
ax.set_ylabel("Fe species [M]")
plt.savefig('enterobactin_decay_linear.png')


