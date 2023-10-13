import numpy as np
from numba import njit
from time import time, strftime, gmtime
from scipy.constants import e, m_e, k as kB, hbar, mu_0

## This first set of routine computes the energy dispersion and defines the Fermi Surface
##
@njit()
def get_E(kx, ky, s_coeff):
    """Compute and return the in-plane dispersion using a TB representation."""
    
    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters
    
    E1 = t1/2*(np.cos(kx*a) + np.cos(ky*b)) # first nearest neighbour hopping energy
    E2 = t2*np.cos(kx*a)*np.cos(ky*b) # second n. n.
    E3 = t3/2*(np.cos(2*kx*a) + np.cos(2*ky*b)) # third n. n.
    E4 = t4/2*(np.cos(2*kx*a)*np.cos(ky*b) + np.cos(kx*a)*np.cos(2*ky*b)) # fourth n. n.
    E5 = t5*np.cos(2*kx*a)*np.cos(2*ky*b) # fifth n. n.
    
    return mu+E1+E2+E3+E4+E5

@njit()
def get_kyF(kxF, s_coeff):
    """At any given kx, compute the corresponding ky on the FS, if it exists."""
    
    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters
    
    ## Solve the equation E=0 for a given kx on the FS
    x = np.cos(kxF*a)
    
    u = t3 + t4*x + t5*(4*x**2-2)
    v = t1/2 + t2*x + t4*(x**2-1/2)
    w = mu + t1/2*x + t3*(x**2-1) - t4/2*x + t5*(1-2*x**2)
    
    y1 = -v/(2*u) + np.sqrt((v/(2*u))**2-w/u)
    y2 = -v/(2*u) - np.sqrt((v/(2*u))**2-w/u)
    
    ky1 = np.arccos(y1)/b
    ky2 = np.arccos(y2)/b
    
    if np.isnan(ky1): return ky2
    else: return ky1 # return np.nan if it does not exist
    
@njit()
def get_kF_lim_1(s_coeff):
    """Compute and return the extreme of the Fermi surface, that is the kx corresponding to ky=π/b."""
    
    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters
    
    ## Solve the equation E=0 for ky=π/b
    u = t3 - t4 + 2*t5
    v = t1/2 - t2 + t4/2
    w = mu - t1/2 + t4/2 - t5
    
    x1 = -v/(2*u) + np.sqrt((v/(2*u))**2-w/u)
    x2 = -v/(2*u) - np.sqrt((v/(2*u))**2-w/u)
    
    kx1 = np.arccos(x1)/a
    kx2 = np.arccos(x2)/a
    
    if (np.isnan(kx1)) and (np.isnan(kx2)):
        return (0., get_kyF(0., s_coeff))
    else: 
        if np.isnan(kx1): return (kx2, np.pi/b)
        else: return (kx1, np.pi/b)
        
@njit()
def get_kF_lim_2(s_coeff):
    """Compute and return the extreme of the Fermi surface, that is the kx corresponding to ky=0."""
    
    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters
    
    ## Solve the equation E=0 for ky=0
    u = t3 + t4 + 2*t5
    v = t1/2 + t2 + t4/2
    w = mu + t1/2 - t4/2 - t5
    
    x1 = -v/(2*u) + np.sqrt((v/(2*u))**2-w/u)
    x2 = -v/(2*u) - np.sqrt((v/(2*u))**2-w/u)
    
    kx1 = np.arccos(x1)/a
    kx2 = np.arccos(x2)/a
    
    if (np.isnan(kx1)) and (np.isnan(kx2)):
        return (np.pi/a, get_kyF(np.pi/a, s_coeff))
    else: 
        if np.isnan(kx1): return (kx2, 0.)
        else: return (kx1, 0.)
        
@njit()
def get_kF_lim(s_coeff):
    """Get the extremes of the Fermi arc in the first quadrant of the BZ."""
    kx_lim_1, ky_lim_1 = get_kF_lim_1(s_coeff)
    kx_lim_2, ky_lim_2 = get_kF_lim_2(s_coeff)
    
    return (kx_lim_1, ky_lim_1), (kx_lim_2, ky_lim_2)

@njit()
def get_kF(s_coeff, N_k, in_AFBZ):
    """Compute and return the Fermi wavevectors evenly spaced along the FS."""
    
    a, b, c = s_coeff[:3] # lattice parameters
    
    (kx_lim_1, ky_lim_1), (kx_lim_2, ky_lim_2) = get_kF_lim(s_coeff) # the extremes of the Fermi arc
    kxF = np.linspace(kx_lim_1, kx_lim_2, N_k)
    kyF = np.empty_like(kxF)
    for i in range(len(kxF)):
        kyF[i] = get_kyF(kxF[i], s_coeff)
    kyF[0], kyF[-1] = ky_lim_1, ky_lim_2
    
    ## Evenly space the k-points along the Fermi arc
    dkF = np.sqrt(np.diff(kxF)**2 + np.diff(kyF)**2) # segment lengths
    FS_l = np.zeros_like(kxF)
    FS_l[1:] = np.cumsum(dkF) # integrate FS path
    FS_grid = np.linspace(0, FS_l.max(), len(kxF)) # regular spaced path
    
    kxF = np.interp(FS_grid, FS_l, kxF)
    kyF = np.empty_like(kxF)
    for i in range(len(kxF)):
        kyF[i] = get_kyF(kxF[i], s_coeff)
    kyF[0], kyF[-1] = ky_lim_1, ky_lim_2
    
    if in_AFBZ: # map inside the AFBZ
        mask_in_AFBZ = np.abs(kyF)<=(np.pi/a-np.abs(kxF))
        kxF, kyF = kxF[mask_in_AFBZ], kyF[mask_in_AFBZ]
        
    ## Negative kxF are obtained by symmetry
    z_ind = 1 if 0. in kxF else 0 # avoid double zero point
    kxF = [-kx for kx in kxF[::-1]] + [kx for kx in kxF[z_ind:]]
    kyF = [ky for ky in kyF[::-1]] + [ky for ky in kyF[z_ind:]]
    
    return np.array(kxF), np.array(kyF)

@njit()
def get_v0(kx, ky, s_coeff):
    """Compute and return the velocities computed as the gradient of the dispersion relation in k-space."""
    
    a, b, c = s_coeff[:3] # lattice constants
    mu, t1, t2, t3, t4, t5 = s_coeff[3:9] # TB parameters
    
    vx = a*np.sin(kx*a)*(-t1-2*t2*np.cos(ky*b)-t4*np.cos(2*ky*b)\
                                            -4*np.cos(kx*a)*(t3+t4*np.cos(ky*b)+2*t5*np.cos(2*ky*b)))/2/hbar
    
    vy = b*np.sin(ky*b)*(-t1-2*t2*np.cos(kx*a)-t4*np.cos(2*kx*a)\
                                                -4*np.cos(ky*b)*(t3+t4*np.cos(kx*a)+2*t5*np.cos(2*kx*a)))/2/hbar
    return vx, vy

@njit()
def get_DOS(kx, ky, s_coeff):
    """Compute and return the density of states."""
    
    vx, vy = get_v0(kx, ky, s_coeff)
    return 1/np.sqrt(vx**2+vy**2)/hbar


## Here, cyclotron motion around the FS is implemented
##
@njit()
def map_inside_BZ(kx, ky, s_coeff, in_AFBZ):
    """Return the wavevectors mapped inside the first BZ (or the AFBZ)."""
    
    a, b, c = s_coeff[:3] # lattice constants
    
    if not in_AFBZ:
        if ky > np.pi/b: ky -= 2*np.pi/b
        if ky < -np.pi/b: ky += 2*np.pi/b
        if kx > np.pi/a: kx -= 2*np.pi/a
        if kx < -np.pi/a: kx += 2*np.pi/a
            
    if in_AFBZ:
        if ky > (np.pi/a+kx): kx, ky = kx+np.pi/a, ky-np.pi/b 
        if ky > (np.pi/a-kx): kx, ky = kx-np.pi/a, ky-np.pi/b
        if ky < (kx-np.pi/a): kx, ky = kx-np.pi/a, ky+np.pi/b
        if ky < (-kx-np.pi/a): kx, ky = kx+np.pi/a, ky+np.pi/b
    
    return kx, ky

@njit()
def evolve_kF(kx_0, ky_0, B, dt, N_t, s_coeff, in_AFBZ):
    """Evolve the Fermi wavevector along the FS according tp cyclotron motion 
    using Runge-Kutta 4-th order algorithm."""
    
    a, b, c = s_coeff[:3] # lattice constants
    
    kx_t, ky_t = np.zeros(N_t), np.zeros(N_t)
    kx_t[0], ky_t[0] = kx_0, ky_0
    
    dk = e*B/hbar
        
    for t in range(1, N_t):
        vx, vy = get_v0(kx_t[t-1], ky_t[t-1], s_coeff)
        K1_x, K1_y = dk*vy, -dk*vx

        vx, vy = get_v0(kx_t[t-1]+dt*K1_x/2, ky_t[t-1]+dt*K1_y/2, s_coeff)
        K2_x, K2_y = dk*vy, -dk*vx

        vx, vy = get_v0(kx_t[t-1]+dt*K2_x/2, ky_t[t-1]+dt*K2_y/2, s_coeff)
        K3_x, K3_y = dk*vy, -dk*vx

        vx, vy = get_v0(kx_t[t-1]+dt*K3_x, ky_t[t-1]+dt*K3_y, s_coeff)
        K4_x, K4_y = dk*vy, -dk*vx

        kx = kx_t[t-1] + dt/6*(K1_x + 2*K2_x + 2*K3_x + K4_x)
        ky = ky_t[t-1] + dt/6*(K1_y + 2*K2_y + 2*K3_y + K4_y)
        
        kx_t[t], ky_t[t] = map_inside_BZ(kx, ky, s_coeff, in_AFBZ) # at every step check if it is inside the BZ
            
    return kx_t, ky_t


## These routines define the scattering rate
##
@njit()
def get_phi(kx, ky, s_coeff, hole_like):
    """Returns the azimuthal angle φ."""
    
    a, b, c = s_coeff[:3] # lattice constants
    
    if hole_like==False: # centered at Γ
        return np.arctan2(ky, kx)
    
    elif hole_like==True: # centered at (π,π)
        x = (np.pi/a-np.abs(kx))
        y = (np.pi/b-np.abs(ky))
        phi = (-np.pi/2-np.arctan(x/y)) * ((kx>=0)&(ky>=0))
        phi += -np.arctan(y/x) * ((kx<0)&(ky>=0))
        phi += np.arctan(y/x) * ((kx<0)&(ky<0))
        phi += (np.pi/2+np.arctan(x/y)) * ((kx>=0)&(ky<0))
        return phi
    
@njit()
def get_tau_inv(kx, ky, T, s_coeff, tau_choice, tau_coeff, hole_like):
    """Compute and return the scattering rate at a given azimuthal angle phi at temperature T."""
    
    if tau_choice[:8] == 'additive': # iso + anisotropic part
        phi = get_phi(kx, ky, s_coeff, hole_like)
        phase = tau_coeff[8]
        
        iso = tau_coeff[0] + tau_coeff[1]*kB*T/hbar + tau_coeff[2]*T**2
        ani = tau_coeff[3] + tau_coeff[4]*T + tau_coeff[5]*T**2
        return iso + ani*np.absolute(np.cos(2*(phi-phase)))**tau_coeff[6]
    
    elif tau_choice == 'fixed peaks': # anisotropic with the peaks intensity is kept fixed
        phi = get_phi(kx, ky, s_coeff, hole_like)
        phase = tau_coeff[8]
        return (tau_coeff[0]*T)/(tau_coeff[7]*(1-np.absolute(np.cos(2*(phi-phase)))**tau_coeff[6])+T)
    
    elif tau_choice == 'DOS': # iso + anisotropic part proportional to the DOS
        dos = get_DOS(kx, ky, s_coeff)
        
        iso = tau_coeff[0] + tau_coeff[1]*kB*T/hbar + tau_coeff[2]*T**2
        ani = tau_coeff[3] + tau_coeff[4]*T + tau_coeff[5]*T**2
        return iso + ani*dos/1e30
    
    else:
        raise ValueError("Not a valid value for `tau_choice`; " +\
                         "supported values are `additive`, `fixed peaks` or `DOS`.")
        
## FLEX-CVC approximation
##
@njit()
def get_alpha(kx, ky, s_coeff):
    """Function needed to implement the FLEX-CVC approximation."""
    
    a, b, c = s_coeff[:3] # lattice constants
    
    dist_AFBZ = np.abs(np.abs(kx)+np.abs(ky)-np.pi/b)/np.sqrt(2) # distance from the AFBZ
    
    sig = np.pi/(a*np.sqrt(2)) # the distance is normalized by half the diagonal of the BZ
    alp = 1-dist_AFBZ/sig
    
    kF_lim_1 = get_kF_lim_1(s_coeff)
    kF_lim_2 = get_kF_lim_2(s_coeff)
    alp = alp*(1 - (kx/kF_lim_2[0])**16) * (1-(ky/kF_lim_1[1])**16) # 0 at the extremes of the FS
    
    return alp*0.9

@njit()
def get_v(kx, ky, s_coeff, FLEX):
    """Compute the velocities."""
    
    vx0, vy0 = get_v0(kx, ky, s_coeff)
    
    if FLEX==False: 
        return vx0, vy0
    
    else:
        a, b, c = s_coeff[:3] # lattice constants
        
        alpha = get_alpha(kx, ky, s_coeff)  # Calculate alpha_k
        
        vxp1, vyp1 = get_v0(-ky, -kx, s_coeff)
        vxp2, vyp2 = get_v0(ky, kx, s_coeff)
        vxp, vyp = vxp1*((kx*ky)>=0)+vxp2*((kx*ky)<0), vyp1*((kx*ky)>=0)+vyp2*((kx*ky)<0) # points at opposite sides of the BZ
        
        Jx, Jy = (vx0 + alpha*vxp)/(1-alpha**2), (vy0 + alpha*vyp)/(1-alpha**2)
        return Jx, Jy
    
## Eventually, this last set of routines computes the conductivity tensor using SCTIF formalism for the BTE. Then to obtain the resistivity, this tensor is inverted.
##
@njit()
def compute_time_integral(kxF, kyF, B, T, j, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, hole_like):
    """Compute the time integral in the Shockley-Chambers formula for the conductivity. The integral
    is evaluated using trapezoidal rule, while the quasi-particle move along the FS in a cyclotron orbit."""
    
    N_t = int(t_lim/dt)+1
    
    kx_t, ky_t = evolve_kF(kxF, kyF, -B, dt, N_t, s_coeff, in_AFBZ)
    vx_t, vy_t = get_v0(kx_t, ky_t, s_coeff)
    vj = vx_t if j == 'x' else vy_t
    
    ## The exponential time integral is computed with the trapezoidal rule simultaneously with the other integral
    ## since the cyclotron motion is the same (the timestep is the same).
    tau_inv = get_tau_inv(kx_t, ky_t, T, s_coeff, tau_choice, tau_coeff, hole_like) # scattering rate
    exp_int = np.zeros(N_t)
    for i in range(1,N_t):
        exp_int[i] = exp_int[i-1]+tau_inv[i-1]+tau_inv[i]
    
    return np.trapz(y=vj*np.exp(-exp_int*dt/2), dx=dt)

@njit()
def compute_sigma_ij(B, T, ij, N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    """Compute a single term of the conductivity tensor using Shockley-Chambers formula. The FS integral
    is evaluated using trapezoidal rule, while the quasi-particle move along the FS in a cyclotron orbit."""
    
    a, b, c = s_coeff[:3] # lattice constants
  
    kxF, kyF = get_kF(s_coeff, N_k, in_AFBZ) # subdivide the FS with evenly spaced points
    vxF, vyF = get_v(kxF, kyF, s_coeff, FLEX)
    
    vi = vxF if ij[0] == 'x' else vyF
    vF = np.sqrt(vxF**2+vyF**2)
    
    ## Time integral
    hole_like = 0. not in kxF # `False` if centered at Γ 
    time_int = np.zeros_like(kxF)
    for i in range(len(kxF)):
        time_int[i] = compute_time_integral(kxF[i], kyF[i], B, T, ij[1], dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, hole_like)
    
    ## FS integral
    integrand = (vi/vF)*time_int
    dkx, dky = np.diff(kxF), np.diff(kyF)
    dkF = np.sqrt(dkx**2+dky**2)
    
    ## Trapezoidal rule for non-uniform grid spacing
    tot = 0
    if not hole_like:
        for i in range(1, len(kxF)): tot += (integrand[i]+integrand[i-1])*dkF[i-1]/2
    elif hole_like:
        for i in range(1, len(kxF[:N_k])): tot += (integrand[:N_k][i]+integrand[:N_k][i-1])*dkF[:N_k][i-1]/2
        for i in range(1, len(kxF[N_k:])): tot += (integrand[N_k:][i]+integrand[N_k:][i-1])*dkF[N_k:][i-1]/2
    
    prefactor = e**2/(2*np.pi**2*c*hbar)
    return prefactor * 2 * tot # a factor of two is added to account also for negative ky

@njit()
def compute_resistivity(B, T, N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    """Compute the resistivity (in μΩcm) by inverting the conductivity tensor."""
    
    ## Only two components are needed, while the others follow from symmetry arguments.
    sigma_xx = compute_sigma_ij(B, T, 'xx', N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    sigma_xy = compute_sigma_ij(B, T, 'xy', N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    det_sigma = sigma_xx**2+sigma_xy**2
    
    rho_xx = sigma_xx/det_sigma/1e-8 # in μΩcm
    rho_xy = -sigma_xy/det_sigma/1e-8 # in μΩcm
    
    return rho_xx, rho_xy, -rho_xy, rho_xx

@njit()
def get_tau_iso(T, s_coeff, tau_choice, tau_coeff, in_AFBZ):
    """Return the minimum value of the scattering time. Useful for estimate the t_lim parameters in the computation."""
    
    kx, ky = get_kF(s_coeff, 11, in_AFBZ)
    hole_like = 0. not in kx
    return np.min(get_tau_inv(kx, ky, T, s_coeff, tau_choice, tau_coeff, hole_like)**(-1))

def compute_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_B = np.zeros((4, len(B), len(TB)))
    for i in range(len(B)):
        for j in range(len(TB)):
            tau = get_tau_iso(TB[j], s_coeff, tau_choice, tau_coeff, in_AFBZ)
            t_lim = tau*N_tau
            rho_B[0,i,j], rho_B[1,i,j], rho_B[2,i,j], rho_B[3,i,j] = \
                compute_resistivity(B[i], TB[j], N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    return rho_B

def compute_error_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_err_xx, rho_err_RH = np.zeros(6), np.zeros(6)
    tau = get_tau_iso(TB[0], s_coeff, tau_choice, tau_coeff, in_AFBZ)
    t_lim = tau*N_tau
    
    # absolute error in FS spacing grid
    rho_err_xx[0], rho_err_RH[0], _, _ = \
    compute_resistivity(B[int(B[0]==0)], TB[0], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    rho_err_xx[3], rho_err_RH[3], _, _ = \
    compute_resistivity(B[-1], TB[0], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    # absolute error in the timestep
    rho_err_xx[1], rho_err_RH[1], _, _ = \
    compute_resistivity(B[int(B[0]==0)], TB[0], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    rho_err_xx[4], rho_err_RH[4], _, _ = \
    compute_resistivity(B[-1], TB[0], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    # absolute error in the truncation of the time integral
    rho_err_xx[2], rho_err_RH[2], _, _ = \
    compute_resistivity(B[int(B[0]==0)], TB[0], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    rho_err_xx[5], rho_err_RH[5], _, _ = \
    compute_resistivity(B[-1], TB[0], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    return rho_err_xx, rho_err_RH
    
def B_sweep(B, TB, N_k, dt, N_tau, out_fileB, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    print('Computing magnetoresistance...')
    t0 = time()
    rho_B = compute_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    tot_t = time()-t0
    print(f'  Elapsed time: {strftime("%Hh:%Mm:%Ss", gmtime(tot_t))}\n')
    
    print('Estimating the errors...')
    rho_err_xx, rho_err_RH = compute_error_B_sweep(B, TB, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    rho_err_xx[0] = abs(rho_err_xx[0]-rho_B[0,int(B[0]==0),0])
    rho_err_xx[1] = abs(rho_err_xx[1]-rho_B[0,int(B[0]==0),0])
    rho_err_xx[2] = abs(rho_err_xx[2]-rho_B[0,int(B[0]==0),0])
    rho_err_xx[3] = abs(rho_err_xx[3]-rho_B[0,-1,0])
    rho_err_xx[4] = abs(rho_err_xx[4]-rho_B[0,-1,0])
    rho_err_xx[5] = abs(rho_err_xx[5]-rho_B[0,-1,0])
    rho_err_RH[0] = abs(rho_err_RH[0]-rho_B[1,int(B[0]==0),0])*10/B[int(B[0]==0)]
    rho_err_RH[1] = abs(rho_err_RH[1]-rho_B[1,int(B[0]==0),0])*10/B[int(B[0]==0)]
    rho_err_RH[2] = abs(rho_err_RH[2]-rho_B[1,int(B[0]==0),0])*10/B[int(B[0]==0)]
    rho_err_RH[3] = abs(rho_err_RH[3]-rho_B[1,-1,0])*10/B[-1]
    rho_err_RH[4] = abs(rho_err_RH[4]-rho_B[1,-1,0])*10/B[-1]
    rho_err_RH[5] = abs(rho_err_RH[5]-rho_B[1,-1,0])*10/B[-1]
    print(f'  ρ_xx absolute error (in μΩcm) at the lowest field at T={TB[0]}K ', end='')
    print(f'- dk: {rho_err_xx[0]:.2E}, dt: {rho_err_xx[1]:.2E}, tlim: {rho_err_xx[2]:.2E}')
    print(f'  ρ_xx absolute error (in μΩcm) at the highest field at T={TB[0]}K ', end='')
    print(f'- dk: {rho_err_xx[3]:.2E}, dt: {rho_err_xx[4]:.2E}, tlim: {rho_err_xx[5]:.2E}\n')
    print(f'  R_H absolute error (in mm^3/C) at the lowest field at T={TB[0]}K ', end='')
    print(f'- dk: {rho_err_RH[0]:.2E}, dt: {rho_err_RH[1]:.2E}, tlim: {rho_err_RH[2]:.2E}')
    print(f'  R_H absolute error (in mm^3/C) at the highest field at T={TB[0]}K ', end='')
    print(f'- dk: {rho_err_RH[3]:.2E}, dt: {rho_err_RH[4]:.2E}, tlim: {rho_err_RH[5]:.2E}\n')

    np.savez(out_fileB, B=B, rho_B=rho_B, T=TB)
    print(f'The output has been written to the file {out_fileB}.npz\n')
    return rho_B

def compute_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_T = np.zeros((4, len(BT), len(T)))
    for j in range(len(T)):
        tau = get_tau_iso(T[j], s_coeff, tau_choice, tau_coeff, in_AFBZ)
        t_lim = tau*N_tau
        for i in range(len(BT)):
            rho_T[0,i,j], rho_T[1,i,j], rho_T[2,i,j], rho_T[3,i,j] = \
                compute_resistivity(BT[i], T[j], N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    return rho_T
    
def compute_error_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    rho_err_xx, rho_err_RH = np.zeros(6), np.zeros(6)
    tau = get_tau_iso(T[0], s_coeff, tau_choice, tau_coeff, in_AFBZ)
    t_lim = tau*N_tau
    
    # absolute error in FS spacing grid
    rho_err_xx[0], rho_err_RH[0], _, _ = \
    compute_resistivity(BT[-1], T[0], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    # absolute error in the timestep
    rho_err_xx[1], rho_err_RH[1], _, _ = \
    compute_resistivity(BT[-1], T[0], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    # absolute error in the truncation of the time integral
    rho_err_xx[2], rho_err_RH[2], _, _ = \
    compute_resistivity(BT[-1], T[0], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    tau = get_tau_iso(T[-1], s_coeff, tau_choice, tau_coeff, in_AFBZ)
    t_lim = tau*N_tau
    
    rho_err_xx[3], rho_err_RH[3], _, _ = \
    compute_resistivity(BT[-1], T[-1], 3*N_k, dt, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    rho_err_xx[4], rho_err_RH[4], _, _ = \
    compute_resistivity(BT[-1], T[-1], N_k, dt/2, t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    
    rho_err_xx[5], rho_err_RH[5], _, _ = \
    compute_resistivity(BT[-1], T[-1], N_k, dt, 5*t_lim, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)

    return rho_err_xx, rho_err_RH
    
def T_sweep(T, BT, N_k, dt, N_tau, out_fileT, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX):
    print('Computing resistivity vs T...')
    t0 = time()
    rho_T = compute_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
    tot_t = time()-t0
    print(f'  Elapsed time: {strftime("%Hh:%Mm:%Ss", gmtime(tot_t))}\n')
    
    if BT[-1] != 0:
        print('Estimating the errors...')
        rho_err_xx, rho_err_RH = compute_error_T_sweep(T, BT, N_k, dt, N_tau, s_coeff, tau_choice, tau_coeff, in_AFBZ, FLEX)
        rho_err_xx[0] = abs(rho_err_xx[0]-rho_T[0,-1,0])
        rho_err_xx[1] = abs(rho_err_xx[1]-rho_T[0,-1,0])
        rho_err_xx[2] = abs(rho_err_xx[2]-rho_T[0,-1,0])
        rho_err_xx[3] = abs(rho_err_xx[3]-rho_T[0,-1,-1])
        rho_err_xx[4] = abs(rho_err_xx[4]-rho_T[0,-1,-1])
        rho_err_xx[5] = abs(rho_err_xx[5]-rho_T[0,-1,-1])
        rho_err_RH[0] = abs(rho_err_RH[0]-rho_T[1,-1,0])*10/BT[-1]
        rho_err_RH[1] = abs(rho_err_RH[1]-rho_T[1,-1,0])*10/BT[-1]
        rho_err_RH[2] = abs(rho_err_RH[2]-rho_T[1,-1,0])*10/BT[-1]
        rho_err_RH[3] = abs(rho_err_RH[3]-rho_T[1,-1,-1])*10/BT[-1]
        rho_err_RH[4] = abs(rho_err_RH[4]-rho_T[1,-1,-1])*10/BT[-1]
        rho_err_RH[5] = abs(rho_err_RH[5]-rho_T[1,-1,-1])*10/BT[-1]
        print(f'  ρ_xx absolute error (in μΩcm) at the lowest temperature at B={BT[-1]}T ', end='')
        print(f'- dk: {rho_err_xx[0]:.2E}, dt: {rho_err_xx[1]:.2E}, tlim: {rho_err_xx[2]:.2E}')
        print(f'  ρ_xx absolute error (in μΩcm) at the highest temperature at B={BT[-1]}T ', end='')
        print(f'- dk: {rho_err_xx[3]:.2E}, dt: {rho_err_xx[4]:.2E}, tlim: {rho_err_xx[5]:.2E}\n')
        print(f'  R_H absolute error (in mm^3/C) at the lowest temperature at B={BT[-1]}T ', end='')
        print(f'- dk: {rho_err_RH[0]:.2E}, dt: {rho_err_RH[1]:.2E}, tlim: {rho_err_RH[2]:.2E}')
        print(f'  R_H absolute error (in mm^3/C) at the highest temperature at B={BT[-1]}T ', end='')
        print(f'- dk: {rho_err_RH[3]:.2E}, dt: {rho_err_RH[4]:.2E}, tlim: {rho_err_RH[5]:.2E}\n')

    np.savez(out_fileT, B=BT, rho_T=rho_T, T=T)
    print(f'The output has been written to the file {out_fileT}\n')
    return rho_T