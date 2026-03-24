## Code accompanying the manuscript by McRose et al., submitted to ISME Communications in 2025

# Requirements for running the model 

The code is pure python and should run with any version of Python
3.X. It makes use of the following Python packages: SciPy, NumPy, and
matplotlib, so those should be installed.

The main code, containg the settings of all reaction rates for the
different Fe-binding ligands, and the definition of the model
equations, is contained in the file `McRoseMorel_ISMEcom_solve.py`.
You can run an example version of the model by just typing

``` python McRoseMorel_ISMEcom_solve.py ```

or (depending on your Python installation) by just clicking on the
file in some graphical user interface. This will lead you through a
few questions and produce one model output plot.

If you want to play around with the model in a more free way, you can
use the model code as a module and write your own calls to the
models. An example for that is done in
'Fig_different_rates_FerrioxiamineB.py' which produces a plot that
compares how the different ligand reaction rates for FerrioxiamineB
from Witter et al, and from Boiteau et al. affect the model outcome.

The code contained in `McRoseMorel_ISMEcom_solve.py` actually contains
four different kinetic models of increasing complexity: 

- The first of these (coded in function 'modelEntB') is for just solving differential equations (1) and (2) for Fe' and for FeY from the manuscript, assuming that the concentration of FeEDTA is constant, and that the concentration of free siderophore Y can be calculated from added siderophore minus FeY. 

- The second (function 'model2EntB') solves the same equation, but adds one more equation describing the decrease of FeEDTA with time, equation (4) in the manuscript.

- The third (function 'model2EntB') extends the second by adding one more equation for the concentration of free ligand, equation (3) from  the manuscript. Adding this equation allows to also describe a degradation of the app form of the added siderophore with time.

- And finally, the fourth adds a description of biological uptake  of Fe' by solving a fifth equation describing the change of bacterial cell numbers.
