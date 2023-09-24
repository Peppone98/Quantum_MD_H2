# Quantum molecular dynamics of the hydrogen molecule

## Run the code

The user has to compile the Makefile by prompting `make` on a terminal window in the same folder of the Makefile. Then the program can be executed with `./H2 OPTION_NAME`. 

The possible choices for `OPTION_NAME` are the following:
* Programs in which Hartree-Fock energy is computed: `SC_HF, BOMD_HF, CPMD_HF, CGMD_HF`.
* Programs in which DFT energy is computed: `SC_DFT, BOMD_DFT, CPMD_DFT, CGMD_DFT`.
* Car-Parrinello evolution of orbital energy: `Evolve_coefficients`.
* Scalar product between sets of coefficients: `Scal_prod`.
* Print the integrals of the exchange and correlation matrix: `EX_CORR`

## Main outputs
The code returns one coefficients file with the name `AAMD_BBB_coeff.txt` and a file with statistics `AAMD_BBB_X_energies.txt`, where `AA` can be `BO`, `CP` or `CG` and `BBB` can be either `HF` or `DFT`.


The columns of the `X_energies.txt` files are organized as follows: `X` (internuclear distance), `E_0` (HF or DFT energy), `T_N` (nuclear kinetic energy), `E_ee` (electron-electron energy), `E_ob` (energy associated with one-body operators), `E_xc` (exchange and correlation energy), `T_f` fictitious kinetic energy, `F` (force given by the gradient of the energy).

```
X          E_0         T_N           E_ee        E_ob        E_xc          T_f            F
1.43462    -1.17521    0.0218615    -0.635305    -2.47389    -0.0985913    1.55495e-06    0.259554
1.45561    -1.17467    0.0213696    -0.631927    -2.45949    -0.0983689    1.48501e-06    0.264643
1.47635    -1.17401    0.0206266    -0.628632    -2.44545    -0.0981512    1.40099e-06    0.269236
1.49673    -1.17324    0.0196661    -0.625436    -2.43185    -0.0979398    1.3063e-06    0.273357
1.51663    -1.17237    0.0185238    -0.622361    -2.41874    -0.0977344    1.20406e-06    0.277013
```

The fictitious kinetic energy `T_f` column is provided only if `CPMD_HF` and `CPMD_DFT` are specified. Moreover, the column containing the exchange and correlation energy `E_xc` is reported only in `DFT` runs.

The user can modify the list of parameters of the program in the header file `definitions.h`:

```
#define FUNCTIONAL_X XC_LDA_X
#define FUNCTIONAL_C XC_LDA_C_PZ

const int iter = 800; /* iterations used in BOMD or Conjugate Gradient */
const double m = 2.0; /* fictitious mass for electronic problem */
const double gamma_el= 1.0; /* electronic damping */
const double M_N = 1836.5; /* nuclear mass */
const double gamma_N = 15.0; /* nuclear damping */
const double h = 0.1; /* electronic time scale */
const double h_N = 43*h; /* nuclear time scale*/
const int N = 4; /* basis centerd on each atom */
const double a[N] = {13.00773, 1.962079, 0.444529, 0.1219492}; /* exponents of the Gaussians */
const double pi = 3.141592653589793; 
const double a_x = 0.0; /* for inclusion of exchange functional */
```
 