# COMPUTING CUPRATES TRANSPORT PROPERTIES BY SOLVING BOLTZMANN TRANSPORT EQUATION
## *I* - INTRODUCTION. THE CUPRATES
### ELECTRONIC STRUCTURE

The structure common to all cuprates is a stack of 2D conducting CuO<sub>2</sub> planes separated by charge reservoir planes. 

<a id='img-cup'><table>
<caption ><h6><i>Fig. 1<i> - Layered copper oxides are composed of CuO<sub>2</sub> planes, typically separated by insulating spacer layers. The electronic structure of these planes primarily involves hybridization of a $3d_{x^2-y^2}$ hole on the copper sites with planar-coordinated $2p_x$ and $2p_y$ oxygen orbitals<a href="#ref-Kei"><sup>[1]</sup></a>.</h6></caption>
<tbody><tr><td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/ff3bc961-41fe-4ad9-b2e1-2430deda148c" width="500"/></td></tr></tbody>
</table></a>

Since the copper atom has a valence of $+2$, its valence shell contains 9 electrons in the $3d$ orbital (or a single hole). A combination of Jahn-Teller distortion and crystal field splitting breaks the $m$-degeneracy of the $d$ orbitals, and the single hole is located in the highest energy orbital, i.e. the $3d_{x^2-y^2}$. This hole is well-localized on the Cu atom due to the strong Coulomb interactions with the $2p$ orbitals of the neighbouring O atoms, so that it takes a large energy to remove an electron from one site and add it to another. This effect produces a "traffic jam" of electrons<a href="#ref-Kei"><sup>[1]</sup></a>. An insulator produced by this classical jamming effect is referred to as a "Mott insulator". However, even a localized electron has a spin whose orientation remains a dynamical degree of freedom. Virtual hopping of these electrons produces, via the Pauli exclusion principle, an antiferromagnetic interaction between neighbouring spins, through the O atom placed between them. This, in turn, leads to a simple antiferromagnetic (Néel) ordered phase below room temperature, in which there are static magnetic moments on the Cu sites with a direction that reverses from one Cu to the next. The resulting 2D energy dispersion can be expressed using a tight-binding representation:

<a id='eq-tb'></a>
$$\varepsilon(k_x,k_y) = \mu+\frac{t_1}{2}\left[\cos(k_xa)+\cos(k_yb)\right]+t_2\cos(k_xa)\cos(k_yb)+\frac{t_3}{2}\left[\cos(2k_xa)+\cos(2k_yb)\right] + \tag{1}$$
$$+ \frac{t_4}{2}\left[\cos(2k_xa)\cos(k_yb)+\cos(k_xa)\cos(2k_yb)\right] + t_5\cos(2k_xa)\cos(2k_yb)$$

where $\mu$ is the chemical potential, $t_1$, $t_2$, $t_3$, $t_4$ and $t_5$ are the hopping parameters up to fifth nearest-neighbour interactions and $a$, $b$ and $c$ are the lattice parameters.

<a id='img-FS'></a>
<table>
<caption >
<h6><i>Fig. 2</i> - The computed tight-binding energy dispersion (in eV) and the Fermi surface (left) and the density of states (right) for La<sub>2-x</sub>Sr<sub>x</sub>CuO<sub>4</sub> (LSCO) at a doping value of $p=0.24$.</h6>
</caption>
<tbody>
<tr>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/7b0f1284-8758-4cc4-937c-aff938916c8f" width="750"/></td>
</tr>
</tbody>
</table>

### PHASE DIAGRAM
The Cu-O planes are ‘doped’ by changing the chemical makeup of interleaved ‘charge-reservoir’ layers so that electrons are removed (hole-doped) or added (electron-doped) to the copper oxide planes. We will consider only hole-doped systems. Hole doping rapidly suppresses the antiferromagnetic order. At a critical doping of $p_{min}$, superconductivity sets in, with a transition temperature that grows to a maximum at $p_{opt}$, then declines for higher dopings and vanishes for $p_{max}$. Materials with $p < p_{opt}$ are referred to as underdoped and those with $p > p_{opt}$ are referred to as overdoped.

<a id='img-phd'></a>
<table>
<caption >
<h6><i>Fig. 3</i> - Temperature versus hole doping level for the copper oxides, indicating where various phases occur<a href="#ref-Kei"><sup>[1]</sup></a>.</h6> 
</caption>
<tbody>
<tr>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/a0d4b46b-ad6d-4268-b4c7-45a9c5060c2e" width="450"/></td>
</tr>
</tbody>
</table>

It is important to recognize that the strong electron repulsions that cause the undoped system to be an insulator (with an energy gap of 2 eV) are still the dominant microscopic interactions, even in optimally doped copper oxide superconductors. This has several general consequences. The resulting electron fluid is ‘highly correlated’, in the sense that for an electron to move through the crystal, other electrons must shift to get out of its way. In contrast, in the Fermi liquid description of simple metals, the quasiparticles (which can be thought of as ‘dressed’ electrons) propagate freely through an effective medium defined by the rest of the electrons.
The failure of the quasiparticle paradigm is most acute in the ‘strange metal’ regime, that is, the ‘normal’ state out of which the pseudogap and the superconducting phases emerge when the temperature is lowered. Nonetheless, in some cases, despite the strong correlations, an emergent Fermi liquid arises at low temperatures. This is especially clear in the overdoped regime.

## *II* - MODELLING THE TRANSPORT PROPERTIES. THE SCTIF

We evaluate the following standard expression for the conductivity within a Boltzmann framework, in the **relaxation time approximation (RTA)**, using the Shockley-Chambers tube-integral formula (SCTIF)<a href="#ref-Sho"><sup>[2,3]</sup></a>,

$$\sigma_{ij} = \frac{e^2}{4\pi^3\hbar}\int_{\text{FS}}d^2\boldsymbol{k} \  \frac{v_i}{|\boldsymbol{v}|}\int^{+\infty}_0dt \ v_j(\boldsymbol{k}(-t)) \ e^{-\int_0^t dt'\ \tau^{-1}(\boldsymbol{k}(-t'))}$$

where $i$ and $j$ corresponds to the crystal directions $a,b,c \equiv x,y,z$ and the integral is along the Fermi Surface (FS). The velocity of the quasiparticle $v_i$ is defined through the dispersion relation, $\hbar v_i = \partial \epsilon/\partial k_i$. $\tau$ is the scattering lifetime. 
This expression is simplified for a cylindrical FS, since $\int d^2k = \frac{2\pi}{c}\int dk$. The final formula is thus

<a id='eq-sigma'></a>
$$\sigma_{ij} = \frac{e^2}{2\pi^2c\hbar}\int_{\text{FS}}d\boldsymbol{k} \ \frac{v_i}{v_F}\int_0^{+\infty}dt \ v_j(\boldsymbol{k}(-t))\ e^{-\int_0^t dt'\ \tau^{-1}(\boldsymbol{k}(-t'))} \tag{2}$$

We assume that the quasiparticle follows standard cyclotron motion, and therefore, we compute the time evolution of $\boldsymbol{k}$ using Newton's equations:

<a id='eq-cycl'></a>
$$\hbar \frac{d\boldsymbol{k}}{dt} = e\boldsymbol{v} \times \boldsymbol{B} \tag{3} $$

We take the magnetic field $\boldsymbol{B}$ oriented along the $z$-direction.

<a id='fig-cycl_orb'></a>
<table>
<caption >
<h6><i>Fig. 4</i> - Simulated cyclotron orbit of the quasiparticle around the Fermi Surface.</h6>
</caption>
<tbody>
<tr>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/4f9efdc0-9db0-4012-8217-5dfe433cc181" width="500"/></td>
</tr>
</tbody>
</table>

We model the scattering time $\tau(\phi)$ with a four-fold azimuthal angle dependence consistent with the crystal symmetry. It consists of two terms: an *isotropic* term (angle-independent) and an *anisotropic* term (angle-dependent). We consider the scattering to depend only on temperature and not on the magnetic field.

<a id='fig-tau_inv'></a>
<table>
<caption >
<h6><i>Fig. 5</i> - Polar logarithmic plot of the scattering rate $\tau^{-1}$ as a function of the azimuthal angle $\phi$.</h6>
</caption>
<tbody>
<tr>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/2a920bc6-9871-473b-8f0a-868f07bb41d1" width="400"/></td>
</tr>
</tbody>
</table>

Eventually, we evaluate the resistivity by inverting the conductivity tensor

<a id='eq-rho'></a>
```math
\rho = \frac{1}{\sigma _{xx}^2+\sigma _{xy}^2} \left( \begin{array}{cc} 
\sigma _{xx} & -\sigma _{xy} \\ 
\sigma _{xy} & \sigma _{xx} \\ 
\end{array} \right) \tag{4}
```

since, for symmetry reason, $\sigma_{xx} = \sigma_{yy}$ and $\sigma_{xy} = -\sigma_{yx}$.

### FLEX-CVC APPROXIMATION
We also consider the Fluctuation Exchange Approximation (FLEX) + Current Vertex Correction (CVC)<a href="#ref-Kon"><sup>[4]</sup></a>. In this approach, the conductivity is computed using the current rather than the quasiparticle velocity. Furthermore, spin fluctuations introduce correlations between wavevectors located on opposite sides of the FS that affect the current. As a result, the current is no longer parallel to the velocity, and thus, no longer perpendicular to the FS. This effect is more enhanced in the proximity of the antiferromagnetic zone boundaries.

<a id='fig-CVC'></a>
<table>
<caption >
<h6><i>Fig. 6</i> - The FLEX-CVC approximation<a href="#ref-Kon"><sup>[4]</sup></a>. The left panel highlights the correlations introduced inside the BZ due to spin fluctuations. The right panel shows how the current is modified by this approximation.</h6>
</caption>
<tbody>
<tr>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/f34eba49-4cb9-4cc0-a3bb-c952cd3fbed9" width="250"/></td>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/34c50c1b-b70b-42c6-ac96-22cf679ee006" width="250"/></td>
</tr>
</tbody>
</table>

The FLEX-CVC approximation is implemented in a Boltzmann framework by substituting the velocity $\boldsymbol{v_i}$ with the current in the expression for the conductivity [(2)](#eq-sigma). The expression for the current is:

<a id='eq-CVC'></a>
$$\boldsymbol{J_k} = \frac{1}{1-\alpha_{\boldsymbol{k}}^2}\left(\boldsymbol{v_k}+\alpha_{\boldsymbol{k}}\boldsymbol{v}_{\boldsymbol{k}'}\right) \tag{5}$$

with $\alpha_{\boldsymbol{k}}<1$ taking the maximum value around the hot spots, and where $(k_x', k_y') = (-k_y, -k_x)$ for $k_xk_y>0$ and $(k_y, k_x)$ for $k_xk_y<0$.

<a id = 'fig-mfp' ></a>
<table>
<caption >
<h6><i>Fig. 6</i> - The circulation of the scattering path length $\boldsymbol{l_k}=\boldsymbol{v_k}\tau_{\boldsymbol{k}}$ with and without FLEX-CVC (bottom panels) around the Fermi surface of three different samples of LSCO (top panels), that determines the Hall coefficient<a href="#ref-Ong"><sup>[5]</sup></a>. Intuitively, we can think of the effect of the FLEX-CVC as to introduce an "effective Fermi surface"<a href="#ref-Kon"><sup>[4]</sup></a> which affects the transport properties of the system considered.</h6>
</caption>
<tbody>
<tr>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/0a3910f7-ebcc-4bd4-b2f2-6a532b7b9e47" width="700"/></td>
</tr>
</tbody>
</table>

## *III* - CODE IMPLEMENTATION
All the routines necessary to run the code are located inside the file `BTE_routines.py` in the `Code` folder.

### DISPERSION AND FERMI SURFACE
This set of functions is responsible for computing the energy dispersion [(1)](#eq-tb) and defining the Fermi Surface. Each of these functions takes as input parameter the array `s_coeff`, which contains the lattice constants and the hopping parameters.

- The function `get_E` computes and returns the energy dispersion at a given point in $k$-space.
```python
def get_E(kx, ky, s_coeff)
```
- For any point $k_x$, the function `get_kyF` returns the corresponing $k_y$ on the FS, if it exists, by solving the equation $\varepsilon(k_x, k_y) = 0$.

```python
def get_kyF(kxF, s_coeff)
```

- The function `get_kF` defines the FS (for positive $k_y$). Each Fermi arc contains `N_k` wavevectors evenly spaced along the FS. The last input is the boolean parameter `in_AFBZ`, which controls whether the wavevectors are restricted to the first Brillouin zone (BZ) or to the antiferromagnetic Brillouin zone (AFBZ).

```python
def get_kF(s_coeff, N_k, in_AFBZ)
```

- The function `get_v0` computes the velocities at a given point in $k$-space. For the tight-binding dispersion considered [(1)](#eq-tb), they are given by

```math
\begin{align}
v_x &= \frac{a}{2\hbar}\sin(k_xa) \{ -t_1-2t_2\cos(k_yb)-t_4\cos(2k_yb)-4\cos(k_xa) [t_3+t_4\cos(k_yb)+2t_5\cos(2k_yb)] \} \\
v_y &= \frac{b}{2\hbar}\sin(k_yb) \{ -t_1-2t_2\cos(k_xa)-t_4\cos(2k_xa)-4\cos(k_yb) [t_3+t_4\cos(k_xa)+2t_5\cos(2k_xa)] \}
\end{align}
```

```python
def get_v0(kx, ky, s_coeff)
```

- The function `get_DOS` computes and returns the density of state, computed as $DOS = 1/|\nabla\varepsilon_k| = 1/\hbar v_F$.

```python
def get_DOS(kx, ky, s_coeff)
```

### CYCLOTRON MOTION
- The function `evolve_kF` evolve the wavevector (`kx_0`, `ky_0`) along the FS by computing `N_t` time steps, with timestep `dt`, at a given value `B` of the magnetic field. To evolve the quasiparticle position according to cyclotron motion [(4)](#eq-cycl), the [Runge-Kutta algorithm](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) at fourth order is implemented. 

```python
def evolve_kF(kx_0, ky_0, B, dt, N_t, s_coeff, in_AFBZ)
```

### SCATTERING RATE

- The function `get_phi` returns the azimuthal angle $\phi$ for a given point in $k$-space. For hole-like Fermi Surfaces, which are centered at $(\pi, \pi)$ rather than at $\Gamma$, it is essential to appropriately set the boolean parameter `hole_like` to ensure accurate angle computations.

```python
def get_phi(kx, ky, s_coeff, hole_like)
```

- The function `get_tau_inv` computes and returns the scattering rate $\tau^{-1}$. It allows different parameterizations according to the input parameter `tau_choice`:

    - `additive` - the scattering rate is decomposed into two additive components: an isotropic, temperature-dependent scattering rate, which comprises a quadratic term and a Planckian linear term, along with an anisotropic term. The anisotropic term is also quadratic in $T$ and has a four-fold azimuthal angle dependence: 
    
    $$\tau^{-1}(\phi,T) = \tau^{-1}_{iso} + \alpha _{iso} \frac{k_B T}{\hbar} + \beta _{iso} T^2 + ( \tau^{-1} _{ani} + \alpha _{ani} T + \beta _{ani} T^2) \times |\cos(2\phi)|^{\nu}$$
    
    - `fixed peaks` - in this model the scattering rate has a fixed value $\tau^{-1}_{iso}$ at the ‘anti-nodal’ regions of the Brillouin zone($\phi = 0°$, $\phi = 90°$, $\phi = 180°$ and $\phi = 270°$). The anisotropy decreases with $T$: 
    
    $$\tau^{-1}(\phi,T) = \tau^{-1} _{iso} \frac{T}{c _{ani} \left( 1 - \left| \cos ( 2\phi ) \right| ^\nu \right) +T}$$
    
    - `DOS` - the same as `additive`, but with the anisotropic term multiplied by the density of states rather than the cosine term: 
    
    $$\tau^{-1}(\phi,T) = \tau^{-1} _{iso} + \alpha _{iso} \frac{k_B T}{\hbar} + \beta _{iso} T^2 + \left( \tau^{-1} _{ani} + \alpha _{ani} T + \beta _{ani} T^2 \right) \times DOS(\phi) \times 1 eVÅ$$

The coefficients for the scattering rate are given in the array `tau_coeff` in the following order: $\tau^{-1} _{iso},\ \alpha _{iso},\ \beta _{iso},\ \tau^{-1} _{ani},\ \alpha _{ani},\ \beta _{ani},\ \nu,\ c _{ani}$ and $\texttt{phase}$.
The parameter `phase` changes the angle-dependence of $\tau^{-1}$ by performing the rotation $\phi \to \phi-\texttt{phase}$. 

```python
def get_tau_inv(kx, ky, T, s_coeff, tau_choice, tau_coeff, hole_like)
```

### FLEX-CVC APPROXIMATION

- To decide whether or not to include the FLEX-CVC approximation in the calculations, one should set the boolean parameter `FLEX` when calling the `get_v` function. If `False`, the velocities are computed using the standard dispersion relation, while, if set to `True`, the function computes and returns the current as given by [(5)](#eq-CVC).

```python
def get_v(kx, ky, s_coeff, FLEX)
```

- The function $\alpha_{\boldsymbol{k}}$ is computed by the routine `get_alpha`: it is set to $0$ at the anti-nodal regions of the BZ and reaches a maximum value $< 0.9$ at points in $k$-space where the FS intersects the AFBZ.

```python
def get_alpha(kx, ky, s_coeff)
```

### RESISTIVITY
To compute the conductivity $\sigma$ using the SCTIF [(2)](#eq-sigma), we need to numerically evaluate three integrals, which will be performed using the [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule). This numerical integration method is implemented in the [numpy function 'trapz'](https://numpy.org/doc/stable/reference/generated/numpy.trapz.html).

- The function `compute_time_integral` computes both the time integral in $dt$ and in $dt'$: provided that the timesteps `dt` and `dt'` are equal, the time evolution of the wavevector along the FS remains the same, allowing us to compute both time integrals simultaneously.

```python
def compute_time_integral(kxF, kyF, B, T, j, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, hole_like)
```

- The function `compute_sigma_ij` computes and returns the $ij$-th element of the conductivity tensor by performing the integral in $dk$ along the FS. Using symmetry arguments, we limit the calculations to the positive $k_y$ sheet, and to account for the other, we apply a factor of $2$. In this routine the trapezoidal rule is manually implemented, since the differential increment $d\boldsymbol{k}$ is not constant.

```python
def compute_sigma_ij(B, T, ij, N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
```

- The function `compute_resistivity` evaluate the resistivity (in μΩcm) by inverting the conductivity tensor.

```python
def compute_resistivity(B, T, N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
```

To facilitate the computations, the two functions `B_sweep` and `T_sweep` execute, respectively, a field sweep and a temperature sweep. Then, they store the results in a numpy '.npz' file. These sweeps can be performed at various temperature (field) values using the arrays `TB` (`BT`). The upper limit of the time integral is computed as a multiple of the minimum scattering time.

```python
def B_sweep(B, TB, N_k, dt, N_tau, out_fileB, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
```
```python
def T_sweep(T, BT, N_k, dt, N_tau, out_fileT, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
```

These functions also estimate the numerical errors introduced during the computation of the integrals for the resistivity (in $\mu\Omega cm$) and for the Hall effect (in $mm^3/C$):
- The integration error in $N_k$, which corresponds to discretization in $k$-space, is determined as the difference from the resistivity value (in $\mu\Omega cm$) calculated using three times the number of $k$-points.
- The integration error in the timestep $dt$ is estimated as the difference from the resistivity value computed using half of the timestep.
- The error in the upper limit of the time integral $N_\tau$ is calculated as the difference from the resistivity value obtained by increasing the upper limit five times.

## *IV* - COMPUTATIONS
To run the actual simulation, three primary components must be specified: the Fermi Surface, determined by the TB parameterization [(1)](#eq-tb), the scattering rate and the parameters of the computation (number of $k$-points, timestep $dt$, etc.). Note that the code performs all the computations using SI units.

The first step is to define the FS. This is accomplished by specifying the array `s_coeff`. To simplify this process, the file `metals.py` includes a dictionary containing pre-defined TB parameterizations for various cuprates (in eV), along with the `Metal` class. Using this class is straightforward: begin by instantiating the object by selecting the metal type and doping level from the dictionary. The necessary array is then stored in the class property `s_coeff`. Here's an example:

```python
LSCO = Metal("LSCO", p="0.23")
s_coeff = LSCO.s_coeff
```

TB parameterizations are provided for Tl<sub>2</sub>Ba<sub>2</sub>CuO<sub>6+δ</sub> (Tl2201)<a href="#ref-Tam"><sup>[6]</sup></a>, La<sub>2-x</sub>Sr<sub>x</sub>CuO<sub>4</sub> (LSCO)<a href="#ref-Lee"><sup>[7]</sup></a>, Nd-LSCO<a href="#ref-Gri"><sup>[8]</sup></a> and Bi<sub>2</sub>Sr<sub>2</sub>CuO<sub>6+δ</sub> (Bi2201)<a href="#ref-He"><sup>[9]</sup></a>. 
The class also provides funtions for plotting various system properties, such as the energy dispersion and the FS (with the method `plot_FS`), the velocities (`plot_vels`) and the density of states (`plot_DOS`).

The scattering rate, as already discussed, is defined by the parameter `tau_choice` and the array `tau_coeff`. Below is reported an example code for conducting a magnetic field sweep. In this specific case, the parameters of the simulation are taken from <a href="#ref-Gri"><sup>[8]</sup></a>.

```python
from Code.BTE_routines import *
from Code.metals import *
from Code.plot_results import *

NdLSCO = Metal('Nd-LSCO', p='0.24')
s_coeff = NdLSCO.s_coeff

B = np.linspace(0, 90, 40) # 40 points of magnetic field between 0 and 90 T
TB = np.array([4, 40, 100]) # each sweep is done at 4, 40, and 100 K

tau_choice = 'additive'
tau_coeff = np.array([8.65e12, 1.2, 0., 63.5e12, 0., 0., 12, 0., 0.])

N_k = 50 # number of k-points of a single Fermi arc
dt = 1e-15 # timestep (in s)
N_tau = 70 # the upper limit in the time integral is N_tau times the scattering time

in_AFBZ = False
FLEX = False

out_fileB = 'rhoB_NdLSCO'

rho_xx, rho_xy, rho_yx, rho_yy = B_sweep(B, TB, N_k, dt, N_tau, out_fileB, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
```

### PLOTTING THE RESULTS
The results can be visualized with the help of the class `Results` found in the file `plot_results.py`. This class takes as the only input parameter the ".npz" output file; various methods can then be called to plot the results. For example, to plot the magnetoresistance corresponding to the example run reported above:

```python
out_NdLSCO = Results('./rhoB_NdLSCO.npz')

fig, ax = plt.subplots(figsize=(6, 3.5), constrained_layout=True)
out_NdLSCO.plot_MR(ax)
out_NdLSCO.make_legend(fig)
```

<a id = 'fig-rhoB' ></a>
<table>
<caption >
<h6><i>Fig. 7</i> - The output produced by the example run. The MR is plotted using the class <code>Results</code>.</h6>
</caption>
<tbody>
<tr>
<td><img src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/97184264-34fb-467b-a00b-e552c98193cd" width="500"/></td>
</tr>
</tbody>
</table>

### GRAPHICAL USER INTERFACE
To set up a field or temperature sweep, it is possible also to make use of a GUI within a [Jupyter Notebook](https://jupyter.org/) (provided that the [ipywidgets](https://ipywidgets.readthedocs.io/en/stable/) package is installed). The usage is straightforward:

```python
from Code.BTE_GUI import *
%matplotlib widget

gui = GUI()
```

<a id = 'fig-rhoB' ></a>
<table>
<caption >
<h6><i>Fig. 8</i> - A video showing how to use the GUI tool to design and run the simulation.</h6>
</caption>
<tbody>
<tr>
<td><video src="https://github.com/acuoghi99/BTE_cuprates/assets/119339344/8285a61e-993b-4e19-8960-61529872d8dd" autoplay loop width="900"/></td>
</tr>
</tbody>
</table>

---

# BIBLIOGRAPHY

1. <a id='ref-Kei'></a> B. Keimer, S. Kivelson, M. Norman *et al.*, From quantum matter to high-temperature superconductivity in copper oxides, [*Nature* **518**, 179–186 (2015)](https://doi.org/10.1038/nature14165).

2. <a id='ref-Sho'></a> W. Shockley, Effect of magnetic fields on conduction - "Tube integrals", [*Phys. Rev.* **79**, 191 (1950)](https://doi.org/10.1103/PhysRev.79.191.2).

3. <a id='ref-Cha'></a> R. G. Chambers, The kinetic formulation of conduction problems, [*Proc. Phys. Soc. A* **65**, 458 (1952)](https://doi.org/10.1088/0370-1298/65/6/114).

4. <a id='ref-Kon'></a> H. Kontani, Anomalous transport phenomena in Fermi liquids with strong magnetic fluctuations, [*Rep. Prog. Phys.* **71**, 026501 (2008)](https://doi.org/10.1088/0034-4885/71/2/026501).

5. <a id='ref-Ong'></a> N. P. Ong, Geometric interpretation of the weak-field Hall conductivity in two-dimensional metals with arbitrary Fermi surface, [*Phys. Rev. B* **43**, 193 (1991)](https://doi.org/10.1103/PhysRevB.43.193).

6. <a id = 'ref-Tam' ></a> C. C. Tam, M. Zhu, J. Ayres *et al.*, Charge density waves and Fermi surface reconstruction in the clean overdoped cuprate superconductor Tl2Ba2CuO6+δ, [*Nat. Commun.* **13**, 570 (2022)](https://doi.org/10.1038/s41467-022-28124-y).

7. <a id='ref-Lee'></a> N. R. Lee-Hone, J. S. Dodge, and D. M. Broun, Disorder and superfluid density in overdoped cuprate superconductors, [*Phys. Rev. B* **96**, 024501 (2017)](https://doi.org/10.1103/PhysRevB.96.024501).

8. <a id='ref-Gri'></a> G. Grissonnanche, Y. Fang, A. Legros *et al.*, Linear-in temperature resistivity from an isotropic Planckian scattering rate, [*Nature* **595**, 667–672 (2021)](https://doi.org/10.1038/s41586-021-03697-8).

9. <a id='ref-He'></a> Y. He et al., Fermi Surface and Pseudogap Evolution in a Cuprate Superconductor, [*Science* **344**, 608-611 (2014)](https://doi.org/10.1126/science.1248221).
