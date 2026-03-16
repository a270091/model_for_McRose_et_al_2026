***Code accompanying the manuscript by McRose et al., submitted to ISME Communications in 2025***

** Requirements for running the model **

The code is pure python and should run with any version of Python
3.X. It makes use of the following Python packages: SciPy, NumPy, and
matplotlib, so those should be installed.

The main code, containg the settings of all reaction rates for the
different Fe-binding ligands, and the definition of the model
equations, is contained in the file 'McRoseMorel_ISMEcom_solve.py'.
You can run an example version of the model by just typing

``` python McRoseMorel_ISMEcom_solve.py ```

or (depending on your Python installation) by just clicking on the
file in some graphical user interface. This will lead you through a
few quations and produce one model output plot.

If you want to play around with the model in a more free way, you can
use the model code as a module and write your own calls to the
models. An example for that is done in
'Fig_different_rates_FerrioxiamineB.py' which produces a plot that
compares how the different ligand reaction rates for FerrioxiamineB
from Witter et al, and from Boiteau et al. affect the model outcome.
