from .BTE_routines import *
from .metals import metals
from .plot_results import *

from os import mkdir
from ipywidgets import Button, Checkbox, Dropdown, FloatRangeSlider, FloatsInput, FloatSlider, FloatText, HBox, interactive_output, IntSlider, IntText, Label, Layout, Output, RadioButtons, Tab, Text, ToggleButtons, VBox

class GUI():
    def __init__(self, root_dir='./Projects/'):
        print('Initializing the routines...', end='')
        
        self.root_dir = root_dir
        
        plt.close('all')
        self.define_system()
        self.define_tau()
        self.define_computation()
        self.make_buttons()
        self.make_layout()
        out_met = interactive_output(self.choose_met, {'met': self.met_w})
        out_dop = interactive_output(self.choose_doping, {'doping': self.doping_w})
        out_latt = interactive_output(self.plot_AFBZ, {'a_latt': self.lattice_a_w, 'b_latt': self.lattice_b_w})
        out_FS = interactive_output(self.plot_FS, {'mu': self.mu_w, 't1': self.t1_w, 't2': self.t2_w, 
                                                   't3': self.t3_w, 't4': self.t4_w, 't5': self.t5_w})
        out_tau_inv = interactive_output(self.plot_tau_inv, {'tau0': self.tau_inv_iso_w, 'a0': self.alpha_iso_w, 
                                                  'b0': self.beta_iso_w, 'tau1': self.tau_inv_ani_w, 
                                                  'a1': self.alpha_ani_w, 'b1': self.beta_ani_w, 
                                                  'nu': self.nu_w, 'phase': self.phase_w})
        out_tau_choice = interactive_output(self.choose_tau_choice, {'t_c': self.tau_choice_w})
        out_B_or_T = interactive_output(self.what_compute, {'B_or_T': self.rho_dep_w})
        
        self.gui_widg = VBox([Text(value=None, 
                                   placeholder='Project name', 
                                   layout=Layout(width='250px', padding='0px 0px 20px 0px')),  
                              
                              Tab([self.w1,  self.w2,  self.w3],  
                                  titles=['SYSTEM', 'SCATTERING RATE', 'COMPUTATION'])])
        
        _ = compute_resistivity(10.0, 4.0, 5, 1e-15, 1e-14, self.get_s_coeff(), 'additive', self.get_tau_coeff(), False, False)
        print('done.')
        
        display(self.gui_widg, out_met, out_dop, out_latt, out_FS, out_tau_inv, out_tau_choice, out_B_or_T)
        
    def define_system(self):
        self.met_w = Dropdown(options=list(metals.keys())+['custom'])
        self.doping_w = Dropdown(options=list(metals[self.met_w.value].keys())[3:], description='doping')
        
        self.lattice_a_w = FloatText(value=1, step=0.01)
        self.lattice_b_w = FloatText(value=1, step=0.01)
        self.lattice_c_w = FloatText(value=1, step=0.01)

        self.mu_w = FloatSlider(description=r'$\mu \ (eV)$', value=1, min=-3, max=3, step=0.0001)
        self.t1_w = FloatSlider(description=r"$t_1 \ (eV)$", value=1, min=-3, max=3, step=0.0001)
        self.t2_w = FloatSlider(description=r"$t_2 \ (eV)$", value=1, min=-3, max=3, step=0.0001)
        self.t3_w = FloatSlider(description=r'$t_3 \ (eV)$', value=1, min=-3, max=3, step=0.0001)
        self.t4_w = FloatSlider(description=r"$t_4 \ (eV)$", value=1, min=-3, max=3, step=0.0001)
        self.t5_w = FloatSlider(description=r"$t_5 \ (eV)$", value=1, min=-3, max=3, step=0.0001)
        
        self.out_fig_FS = Output()
        with self.out_fig_FS:
            self.fig_FS, self.ax_FS = plt.subplots(figsize=(3.25,3.25), constrained_layout=True)
            self.ax_FS.set_aspect('equal')
            self.fig_FS.canvas.header_visible = False
            self.fig_FS.canvas.footer_visible = False
            self.fig_FS.canvas.toolbar_visible = False
            self.fig_FS.canvas.resizable = False
            tmp_2D = np.array([[-1,1],[-1,1]])
            self.im = self.ax_FS.contourf(tmp_2D, tmp_2D, tmp_2D)
            self.FS = self.ax_FS.contour(tmp_2D, tmp_2D, tmp_2D)
            self.cbar = plt.colorbar(self.im, ax=self.ax_FS, orientation='horizontal', format="%.2f", 
                                     location='top')
            plt.show()
            
        eq = r'$$\begin{align} \varepsilon(k_x,k_y) = '
        eq += r'\mu + \frac{t_1}{2}\left[\cos(k_xa)+\cos(k_yb)\right] &+ t_2\cos(k_xa)\cos(k_yb) + '
        eq += r'\frac{t_3}{2}\left[\cos(2k_xa)+\cos(2k_yb)\right] +\\'
        eq += r'& + \frac{t_4}{2}\left[\cos(2k_xa)\cos(k_yb)+\cos(k_xa)\cos(2k_yb)\right] + '
        eq += r't_5\cos(2k_xa)\cos(2k_yb)\end{align}$$'
        self.eq_tb_w = Label(eq)

    def define_tau(self):
        self.tau_choice_w = Dropdown(description=r'$\tau$ model', 
                                     options=['additive', 'additive (Gr)', 'fixed peaks', 'DOS'])
        self.phase_w = IntSlider(description='phase (°)', value=0, min=0, max=90)
        self.tau_inv_iso_w = FloatSlider(value=10., min=0, max=1e6, step=0.01)
        self.alpha_iso_w = FloatSlider(value=1.0, min=0, max=1e6, step=0.01)
        self.beta_iso_w = FloatSlider(value=0., min=0, max=1e6, step=0.01)
        self.tau_inv_ani_w = FloatSlider(value=0., min=0, max=1e6, step=0.01)
        self.alpha_ani_w = FloatSlider(value=0., min=0, max=1e6, step=0.01)
        self.beta_ani_w = FloatSlider(value=0., min=0, max=1e6, step=0.01)
        self.nu_w = IntSlider(value=2, min=0, max=1000, step=1)
        
        self.out_fig_tau = Output()
        with self.out_fig_tau:
            self.fig_tau, self.ax_tau = plt.subplots(figsize=(3.4,3.), constrained_layout=True)
            self.fig_tau.canvas.header_visible = False
            self.fig_tau.canvas.footer_visible = False
            self.fig_tau.canvas.toolbar_visible = False
            self.fig_tau.canvas.resizable = False
            plt.show()
            
        eq = r'$eq_tau$'
        self.eq_tau_w = Label(eq)
    
    def define_computation(self):
        self.rho_dep_w = ToggleButtons(options=['Field sweep ','Temperature sweep ', 'Both '],
                                       icons=['magnet', 'thermometer-three-quarters', 'expand'], 
                                       button_style='info')
        self.B_w = FloatRangeSlider(description='B range (T)', value=[0., 90], min=0., max=150.0, step=0.1, 
                                    readout_format='.1f')
        self.NB_w = IntText(description='N_B', value=40)
        self.TB_w = FloatsInput(value=[4, 40, 100], readout_format='.1f')
        self.T_w = FloatRangeSlider(description='T range (K)', value=[4, 100], 
                                    min=1, max=150.0, step=0.1, disabled=True, readout_format='.1f')
        self.NT_w = IntText(description='N_T', value=40, disabled=True)
        self.BT_w = FloatsInput(value=[4, 35], disabled=True, readout_format='.1f')
        
        self.Nk_w = IntText(description='N_k', value=25)
        self.dt_w = FloatText(description='dt (fs)', value=1.0, step=0.1)
        self.Ntau_w = IntText(description='N_tau', value=10)
        
        self.inAFBZ_w = Checkbox(value=False, description='inside AFBZ', indent=False)
        self.FLEX_w = Checkbox(value=False, description='FLEX-CVC', indent=False)
        
        self.out_file_w = Text(placeholder='Output file (".npz" file)')
        
    def make_buttons(self):
        self.button_opt_param = Button(description='Optimize', 
                                       tooltip='Optimize the parameters of the computation in order to\n'+\
                                               'find a proper convergence of the integrals. - STILL WORKING ON IT')
        self.button_compute = Button(description='COMPUTE', button_style='success', icon='wrench')
        self.button_compute.on_click(self.compute)
        
    def make_layout(self):
        self.met_w.layout = Layout(width='150px', padding='0px 0px 20px 20px')
        self.doping_w.layout = Layout(width='200px')
        
        lattice_w = HBox([HBox([Label(r'$a \ (Å)$'), self.lattice_a_w]), 
                          HBox([Label(r'$b \ (Å)$'), self.lattice_b_w]), 
                          HBox([Label(r'$c \ (Å)$'), self.lattice_c_w])],
                         layout=Layout(justify_content='space-between', padding='20px 0px 20px 20px'))
        
        for widg in lattice_w.children:
            widg.layout = Layout(justify_content='space-between', align_items='baseline')
            widg.children[1].readout_format = '.2f'
            widg.children[1].layout = Layout(width='100px', padding='0px 20px 0px 10px')
            
        tght_bnd_w = HBox([self.mu_w, self.t1_w, self.t2_w, self.t3_w, self.t4_w, self.t5_w])
        for widg in tght_bnd_w.children:
            widg.readout_format = '.4f'
            widg.orientation = 'vertical'
            widg.layout = Layout(width='90px')
        self.eq_tb_w.layout = Layout(height='100px', padding='5px 0px 0px 100px')

        self.tau_choice_w.layout = Layout(width='200px')
        self.phase_w.layout = Layout(width='370px', padding='0px 0px 0px 70px')
        tau_iso_w = VBox([HBox([Label(r'$\tau^{-1}_{iso}\ (ps^{-1})$'), self.tau_inv_iso_w]), 
                          HBox([Label(r'$\alpha_{iso}$'), self.alpha_iso_w]), 
                          HBox([Label(r'$\beta_{iso}\ (ns^{-1}/K^2)$'), self.beta_iso_w])])
        tau_ani_w = VBox([HBox([Label(r'$\tau^{-1}_{ani}\ (ps^{-1})$'), self.tau_inv_ani_w]),  
                          HBox([Label(r'$\alpha_{ani}\ (ps^{-1}/K)$'), self.alpha_ani_w]), 
                          HBox([Label(r'$\beta_{ani}\ (ns^{-1}/K^2)$'), self.beta_ani_w]), 
                          HBox([Label(r'$\nu$'), self.nu_w])])
        self.tau_inv_w = HBox([tau_iso_w, tau_ani_w], 
                              layout=Layout(justify_content='space-between', padding='40px 0px 0px 0px'))
        
        for widg in self.tau_inv_w.children:
            for wdg in widg.children:
                wdg.layout = Layout(justify_content='flex-end', align_items='baseline')
                wdg.children[1].readout_format = '.2e'
                wdg.children[0].layout = Layout(padding='0px 10px 0px 0px')
                wdg.children[1].layout = Layout(height='60px', width='180px', padding='0px 0px 0px 0px')
        self.nu_w.readout_format = 'd'
        self.eq_tau_w.layout = Layout(height='50px', padding='0px 0px 0px 250px')
        
        self.rho_dep_w.layout = Layout(padding='0px 0px 20px 0px')
        self.rho_dep_w.style.button_width = '200px'
        self.B_w.layout = Layout(width='400px')
        self.NB_w.layout = Layout(width='150px')
        self.T_w.layout = Layout(width='400px')
        self.NT_w.layout = Layout(width='150px')
        
        field_range_w = VBox([HBox([self.B_w, self.NB_w], layout=Layout(padding='10px 0px 0px 21px')),
                              HBox([Label('T (K)'), self.TB_w], layout=Layout(padding='0px 0px 0px 55px')),
                              HBox([self.T_w, self.NT_w], layout=Layout(padding='10px 0px 0px 21px')),
                              HBox([Label('B (T)'), self.BT_w], layout=Layout(padding='0px 0px 0px 55px'))])
        N_comp_w = VBox([self.Nk_w, self.dt_w, self.Ntau_w])
        for widg in N_comp_w.children:
            widg.layout = Layout(width='200px', padding='0px 0px 0px 40px')
        self.Ntau_w.layout = Layout(width='200px', padding='0px 0px 20px 40px')
        self.button_opt_param.layout = Layout(width='200px')
        comp_param_w = VBox([Label('Computation parameters: '), N_comp_w, self.button_opt_param], 
                            layout=Layout(padding='10px 0px 0px 80px', width='400px'))
        self.inAFBZ_w.layout = Layout(width='150px', padding='0px 0px 50px 0px')
        self.FLEX_w.layout = Layout(width='150px', padding='0px 0px 50px 0px')
        self.out_file_w.layout = Layout(padding='0px 0px 10px 0px')
        
        
        self.w1 = VBox([HBox([self.out_fig_FS, VBox([HBox([self.met_w, self.doping_w]), 
                                                     lattice_w, tght_bnd_w])]), self.eq_tb_w])

        self.w2 = VBox([HBox([self.out_fig_tau, VBox([HBox([self.tau_choice_w, self.phase_w]), self.tau_inv_w], 
                             layout=Layout(padding='0px 0px 0px 20px'))]),
                        self.eq_tau_w])

        self.w3 = VBox([self.rho_dep_w, HBox([field_range_w, comp_param_w]), 
                        HBox([self.inAFBZ_w, self.FLEX_w]),
                        self.out_file_w, self.button_compute])
        
    def get_s_coeff(self):
        a, b, c = self.lattice_a_w.value/1e10, self.lattice_b_w.value/1e10, self.lattice_c_w.value/1e10
        mu, t1, t2 = self.mu_w.value*e, self.t1_w.value*e, self.t2_w.value*e
        t3, t4, t5 = self.t3_w.value*e, self.t4_w.value*e, self.t5_w.value*e
        return np.array([a, b, c, mu, t1, t2, t3, t4, t5])
    
    def get_tau_coeff(self):
        tau_inv_iso, alpha_iso, beta_iso = self.tau_inv_iso_w.value, self.alpha_iso_w.value, self.beta_iso_w.value
        tau_inv_ani, alpha_ani, beta_ani = self.tau_inv_ani_w.value, self.alpha_ani_w.value, self.beta_ani_w.value
        nu, phase = self.nu_w.value, self.phase_w.value
        return np.array([tau_inv_iso*1e12, alpha_iso, beta_iso*1e9, 
                         tau_inv_ani*1e12, alpha_ani*1e12, beta_ani*1e9, 
                         nu, tau_inv_ani, phase*np.pi/180])
    
    def choose_met(self, met):
        if met == 'custom':
            self.doping_w.disabled = True
            self.doping_w.value = None
            self.lattice_a_w.disabled, self.lattice_b_w.disabled, self.lattice_c_w.disabled = False, False, False
            self.mu_w.disabled, self.t1_w.disabled, self.t2_w.disabled = False, False, False
            self.t3_w.disabled, self.t4_w.disabled, self.t5_w.disabled = False, False, False
        else:
            self.doping_w.value = None
            self.doping_w.options = list(metals[met].keys())[3:]
            self.doping_w.value = self.doping_w.options[0]
            self.doping_w.disabled = False
            
            self.lattice_a_w.disabled, self.lattice_b_w.disabled, self.lattice_c_w.disabled = True, True, True
            self.mu_w.disabled, self.t1_w.disabled, self.t2_w.disabled = True, True, True
            self.t3_w.disabled, self.t4_w.disabled, self.t5_w.disabled = True, True, True
        
    def choose_doping(self, doping):
        met = self.met_w.value
        if doping != None:
            a, b, c = metals[met]['a'], metals[met]['b'], metals[met]['c']
            mu, t1, t2, t3, t4, t5 = metals[met][doping].values()
            self.lattice_a_w.value, self.lattice_b_w.value, self.lattice_c_w.value = a, b, c
            self.mu_w.value, self.t1_w.value, self.t2_w.value = mu, t1, t2
            self.t3_w.value, self.t4_w.value, self.t5_w.value = t3, t4, t5
        
    def plot_AFBZ(self, a_latt, b_latt):
        self.ax_FS.clear()
        a, b = a_latt*1e-10, b_latt*1e-10
        
        # The antiferromagnetic Brillouin Zone (AFBZ) is half of the standard BZ
        self.ax_FS.plot([-np.pi/a, 0, np.pi/a, 0, -np.pi/a], [0, -np.pi/b, 0, np.pi/b, 0], 
                        '--', color='black', lw=0.6)

        # Plot adjustment
        self.ax_FS.set_xlabel(r'$k_x$', fontsize=11)
        self.ax_FS.set_ylabel(r'$k_y$', fontsize=11)
        self.ax_FS.set(xlim=(-np.pi/a, np.pi/a), ylim=(-np.pi/b, np.pi/b))
        self.ax_FS.set_xticks([-np.pi/a, 0, np.pi/a], labels=[r'$-\pi/a$', 0, r'$\pi/a$'], fontsize=11)
        self.ax_FS.set_yticks([-np.pi/b, 0, np.pi/b], labels=[r'$-\pi/b$', 0, r'$\pi/b$'], fontsize=11)

        self.ax_FS.text(np.pi/a*0.55, np.pi/b*0.55, 'AFBZ', fontsize=10)
        self.ax_FS.plot(0, 0, '.', ms=2., c='black')
        self.ax_FS.text(0, 0, r' $\Gamma$', horizontalalignment='left', verticalalignment='bottom', fontsize=10)
        
        self.mu_w.value += self.mu_w.value*1e-10
        
    def plot_FS(self, mu, t1, t2, t3, t4, t5):
        a, b, c = self.lattice_a_w.value/1e10, self.lattice_b_w.value/1e10, self.lattice_c_w.value/1e10
        s_coeff = np.array([a, b, c, mu*e, t1*e, t2*e, t3*e, t4*e, t5*e])
        kx, ky = np.linspace(-np.pi/a, np.pi/a, 101), np.linspace(-np.pi/b, np.pi/b, 101)
        kx, ky = np.meshgrid(kx, ky)
        energy = get_E(kx, ky, s_coeff)/e # compute the dispersion in eV

        try: 
            for im_coll in self.im.collections: im_coll.remove()
        except:
            pass
        
        try:
            for FS_coll in self.FS.collections: FS_coll.remove()
        except:
            pass
        
        try:
            self.cbar.remove()
        except: 
            pass
        
        self.im = self.ax_FS.contourf(kx, ky, energy, levels=70, cmap='summer')
        self.FS = self.ax_FS.contour(kx, ky, energy, levels=[0], colors=['black'], linewidths=1) # FS

        cbar_ticks = list(np.linspace(np.min(energy[np.isfinite(energy)]), 
                                      np.max(energy[np.isfinite(energy)]), 5))
        self.cbar = plt.colorbar(self.im, ax=self.ax_FS, orientation='horizontal', location='top', 
                                 fraction=0.07, pad=0.1)
        self.cbar.set_label(label=r'$E\ (eV)$', size=10)
        self.cbar.set_ticks(cbar_ticks, labels=[f'{tcks:.2f}' for tcks in cbar_ticks], fontsize=10)
        
        self.tau_inv_iso_w.value += self.tau_inv_iso_w.value*1e-10
            
    def choose_tau_choice(self, t_c):
        if t_c[:8] == 'additive':
            self.tau_inv_w.children[0].children[0].children[0].value = r'$\tau^{-1}_{iso}\ (ps^{-1})$'
            self.tau_inv_w.children[1].children[0].children[0].value = r'$\tau^{-1}_{ani}\ (ps^{-1})$'
            
            eq_add = r'$$\tau^{-1}(\phi,T) = \tau^{-1}_{iso} + \alpha_{iso}\frac{k_B T}{\hbar} +' +\
                     r'\beta_{iso} T^2 + \left(\tau^{-1}_{ani} + \alpha_{ani} T + \beta_{ani} T^2\right)' +\
                     r' \left|\cos(2\phi)\right|^\nu$$'
            self.eq_tau_w.value = eq_add
            
            self.phase_w.disabled = False
            
            if t_c[-1] == ')':
                self.tau_inv_iso_w.disabled = True
                self.alpha_iso_w.disabled = True
                self.beta_iso_w.disabled = True
                self.alpha_ani_w.disabled = True
                self.beta_ani_w.disabled = True
                self.tau_inv_ani_w.disabled = True
                self.nu_w.disabled = True
                
                if t_c[-3:-1] == 'Gr':
                    self.tau_inv_iso_w.value = 8.65
                    self.alpha_iso_w.value = 1.2
                    self.beta_iso_w.value = 0
                    self.tau_inv_ani_w.value = 63.5
                    self.alpha_ani_w.value = 0
                    self.beta_ani_w.value = 0
                    self.nu_w.value = 12
            
            else:
                self.tau_inv_iso_w.disabled = False
                self.alpha_iso_w.disabled = False
                self.beta_iso_w.disabled = False
                self.tau_inv_ani_w.disabled = False
                self.alpha_ani_w.disabled = False
                self.beta_ani_w.disabled = False
                self.nu_w.disabled = False
            
        elif t_c == 'fixed peaks':
            eq_fixd = r'$$\phantom{aaaaaaaaaaaa}\tau^{-1}(\phi,T) = \tau^{-1}_p\frac{T}{c_{ani}'+\
                      r'\left(1-\left|\cos(2\phi)\right|^\nu\right)+T}$$'
            self.eq_tau_w.value = eq_fixd
            
            self.tau_inv_w.children[0].children[0].children[0].value = r'$\tau^{-1}_{p}\ (ps^{-1})$'
            self.tau_inv_w.children[1].children[0].children[0].value = r'$c_{ani}$'
            self.tau_inv_ani_w.value = 7.5
            self.tau_inv_iso_w.disabled = False
            self.alpha_iso_w.disabled = True
            self.beta_iso_w.disabled = True
            self.tau_inv_ani_w.disabled = False
            self.alpha_ani_w.disabled = True
            self.beta_ani_w.disabled = True
            self.nu_w.disabled = False
            self.phase_w.disabled = False
            
        elif t_c == 'DOS':
            eq_DOS = r'$$\tau^{-1}(\phi,T)=\tau^{-1}_{iso}+\alpha_{iso}\frac{k_B T}{\hbar}+\beta_{iso} T^2+'+\
                     r'\left(\tau^{-1}_{ani} + \alpha_{ani} T + \beta_{ani} T^2\right) DOS*1 eVÅ$$'
            self.eq_tau_w.value = eq_DOS
            
            self.tau_inv_w.children[0].children[0].children[0].value = r'$\tau^{-1}_{iso}\ (ps^{-1})$'
            self.tau_inv_w.children[1].children[0].children[0].value = r'$\tau^{-1}_{ani}\ (ps^{-1})$'
            
            self.tau_inv_iso_w.disabled = False
            self.alpha_iso_w.disabled = False
            self.beta_iso_w.disabled = False
            self.alpha_ani_w.disabled = False
            self.beta_ani_w.disabled = False
            self.tau_inv_ani_w.disabled = False
            self.nu_w.disabled = True
            self.phase_w.disabled = True
        self.tau_inv_iso_w.value += self.tau_inv_iso_w.value*1e-10
            
    def plot_tau_inv(self, tau0, a0, b0, tau1, a1, b1, nu, phase):
        cmap = plt.cm.Blues
        color = [cmap(i) for i in np.linspace(0, cmap.N, 5, dtype=np.int32)][2:]
        
        self.ax_tau.clear()
        s_coeff = self.get_s_coeff()
        tau_choice = self.tau_choice_w.value
        in_AFBZ = self.inAFBZ_w.value
        
        tau_coeff = np.array([tau0*1e12, a0, b0*1e9, tau1*1e12, a1*1e12, b1*1e9, nu, tau1, phase*np.pi/180])
        
        N_k = 201
        kx, ky = get_kF(s_coeff, N_k, in_AFBZ)
        hole_like = 0. not in kx
        phi = get_phi(kx, ky, s_coeff, hole_like)
            
        T = [4.0, 40.0, 100.0]
        
        for i in range(3):
            tau_inv = get_tau_inv(kx, ky, T[i], s_coeff, tau_choice, tau_coeff, hole_like)/1e12
            self.ax_tau.plot(phi[:N_k], tau_inv[:N_k], '-', ms=1., c=color[i], label=f'{int(T[i])} K')
            self.ax_tau.plot(phi[N_k:], tau_inv[N_k:], '-', ms=1., c=color[i])
            if phi[0] <= 0:
                self.ax_tau.plot(phi[:N_k]+np.pi, tau_inv[:N_k], '-', ms=1., c=color[i])
                self.ax_tau.plot(phi[N_k:]+np.pi, tau_inv[N_k:], '-', ms=1., c=color[i])
            else:
                self.ax_tau.plot(phi[N_k:]-np.pi, tau_inv[N_k:], '-', ms=1., c=color[i])
                self.ax_tau.plot(phi[:N_k]-np.pi, tau_inv[:N_k], '-', ms=1., c=color[i])

        self.ax_tau.set_xlabel(r'$\phi$', fontsize=11)
        self.ax_tau.set_ylabel(r'$\tau^{-1}\ (ps^{-1})$', fontsize=11)
        self.ax_tau.set(xlim=(-np.pi, np.pi), title=' \n ')
        self.ax_tau.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], 
                               labels=[r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])
        self.ax_tau.tick_params(labelsize=11)
        self.fig_tau.legend(ncol=3, loc='center', bbox_to_anchor=(0.57, 0.95), 
                            fancybox=True, shadow=True, fontsize=10)
        
    def what_compute(self, B_or_T):
        if B_or_T == 'Field sweep ':
            self.B_w.disabled, self.NB_w.disabled, self.TB_w.value = False, False, [4, 40, 100]
            self.T_w.disabled, self.NT_w.disabled, self.BT_w.value = True, True, []
            
        elif B_or_T == 'Temperature sweep ':
            self.B_w.disabled, self.NB_w.disabled, self.TB_w.value = True, True, []
            self.T_w.disabled, self.NT_w.disabled, self.BT_w.value = False, False, [4, 35]
            
        else:
            self.B_w.disabled, self.NB_w.disabled, self.TB_w.value = False, False, [4, 40, 100]
            self.T_w.disabled, self.NT_w.disabled, self.BT_w.value = False, False, [4, 35]
     
    def initialize_dirs(self):
        proj = self.gui_widg.children[0].value
        proj = '.TMP' if proj == '' else proj
        self.proj_dir = self.root_dir+proj+'/'
        self.out_dir = self.proj_dir + 'Output/'
        self.fig_dir = self.proj_dir + 'Figures/'
        for dirs in [self.root_dir, self.proj_dir, self.out_dir, self.fig_dir]:
            try: mkdir(dirs)
            except: pass
     
    def get_params(self):
        self.met = self.met_w.value
        self.doping = self.doping_w.value
        self.s_coeff = self.get_s_coeff()
        self.tau_coeff, self.tau_choice = self.get_tau_coeff(), self.tau_choice_w.value
        self.B_or_T = self.rho_dep_w.value
        self.B_r, self.NB, self.TB = self.B_w.value, self.NB_w.value, sorted(self.TB_w.value)
        self.T_r, self.NT, self.BT = self.T_w.value, self.NT_w.value, sorted(self.BT_w.value)
        self.N_k, self.dt, self.N_tau = self.Nk_w.value, self.dt_w.value/1e15, self.Ntau_w.value
        self.in_AFBZ, self.FLEX = self.inAFBZ_w.value, self.FLEX_w.value
        filename =  self.out_file_w.value
        self.filename = 'rho_TMP' if filename == '' else filename
    
    def print_log(self):
        out_log = f'  - system: {self.met} ({self.doping})\n'
        if self.met == 'custom':
            out_log += f'      * lattice parameters: ' +\
                       f'a = {self.s_coeff[0]:.2E}\n' +\
              ' '*28 + f'b = {self.s_coeff[1]:.2E}\n' +\
              ' '*28 + f'c = {self.s_coeff[2]:.2E}\n' +\
                       f'      * TB parameters: ' +\
                       f'μ = {self.s_coeff[3]:.2E}\n' +\
              ' '*23 + f't1 = {self.s_coeff[4]:.2E}\n' +\
              ' '*23 + f't2 = {self.s_coeff[5]:.2E}\n' +\
              ' '*23 + f't3 = {self.s_coeff[6]:.2E}\n' +\
              ' '*23 + f't4 = {self.s_coeff[7]:.2E}\n' +\
              ' '*23 + f't5 = {self.s_coeff[8]:.2E}\n' 
        out_log += '  - scattering time:\n' +\
                  f'      * tau_choice: {self.tau_choice}\n' +\
                   '      * tau_coeff: '
        if self.tau_choice == 'fixed peaks':
            out_log += f'tau_p = {self.tau_coeff[0]:.2E}\n' +\
              ' '*19 + f'c_ani = {self.tau_coeff[7]:.2f}\n'
        else: 
            out_log += f'tau_iso = {self.tau_coeff[0]:.2E}\n' +\
              ' '*19 + f'alpha_iso = {self.tau_coeff[1]:.2f}\n' +\
              ' '*19 + f'beta_iso = {self.tau_coeff[2]:.2E}\n' +\
              ' '*19 + f'tau_ani = {self.tau_coeff[3]:.2E}\n' +\
              ' '*19 + f'alpha_ani = {self.tau_coeff[4]:.2E}\n' +\
              ' '*19 + f'beta_ani = {self.tau_coeff[5]:.2E}\n'
        if self.tau_choice != 'DOS':
            out_log += ' '*19 + f'nu = {int(self.tau_coeff[6])}\n' +\
                       ' '*19 + f'offset = {int(self.tau_coeff[8]*180/np.pi)}\n'
        out_log += f'  - computation:\n' +\
                   f'      * sweep: {self.B_or_T}\n'
        if (self.B_or_T[0] == 'F') or (self.B_or_T[0] == 'B'):
            out_log += f'          ^ B range: {self.B_r}\n' +\
                       f'          ^ N_B: {self.NB}\n' +\
                       f'          ^ T_B: {self.TB}\n'
        if self.B_or_T[0] == 'B':
            out_log += '\n'
        if (self.B_or_T[0] == 'T') or (self.B_or_T[0] == 'B'):
            out_log += f'          ^ T range: {self.T_r}\n' +\
                       f'          ^ N_T: {self.NT}\n' +\
                       f'          ^ B_T: {self.BT}\n'
        out_log += f'      * N_k: {self.N_k}\n' +\
                   f'      * dt: {self.dt:.2E}\n' +\
                   f'      * N_tau: {self.N_tau}\n' +\
                   f'      * in_AFBZ: {self.in_AFBZ}\n' +\
                   f'      * FLEX-CVC: {self.FLEX}\n' +\
                   f'      * out_file: {self.filename}\n'
        
        with open(self.proj_dir+'output_log.txt', 'a') as f:
            print(f'Run on {strftime("%d-%m-%Y %H:%M:%S GMT", gmtime())}\n', file=f)
            print(out_log, file=f)
            print('#'*50+'\n', file=f)
            
        print(out_log)
        
    def plot_results(self):
        if (self.B_or_T[0] == 'F') or (self.B_or_T[0] == 'B'):   
            Output = Results(self.out_fileB)
            fig_rho = plt.figure(figsize=(9,4.), constrained_layout=True)
            fig_rho.canvas.header_visible = False
            
            ax_rho = fig_rho.add_subplot(2,3,(1,4))
            Output.plot_rho_ij_vs_B('xx', ax_rho, legend=True)
            
            ax_rho = fig_rho.add_subplot(2,3,2)
            Output.plot_MR(ax_rho, xlabel='no', legend=False)
            
            ax_rho = fig_rho.add_subplot(2,3,5)
            Output.plot_RH_vs_B(ax_rho, legend=False)
            
            ax_rho = fig_rho.add_subplot(2,3,3)
            Output.plot_deriv_scaling(ax_rho, legend=False, xlabel='no')
            
            ax_rho = fig_rho.add_subplot(2,3,6)
            Output.plot_Drho_ov_T_scaling(ax_rho, legend=False)
            
            Output.make_legend(fig_rho)
            
            fig_rho.savefig(self.fig_fileB)
            
        if (self.B_or_T[0] == 'T') or (self.B_or_T[0] == 'B'):
            Output = Results(self.out_fileT)
            fig_rho, ax_rho = plt.subplots(1,2, figsize=(8,3), constrained_layout=True)
            fig_rho.canvas.header_visible = False
            
            Output.plot_rho_ij_vs_T('xx', ax_rho[0], legend=True)
            Output.plot_RH_vs_T(ax_rho[1], legend=False, yticks='right')
            Output.make_legend(fig_rho)
            
            fig_rho.savefig(self.fig_fileT)
            
        plt.show()
     
    def compute(self, button_compute):
        print('Starting computation...\n')
        self.initialize_dirs()
        self.get_params()
        self.print_log()
        
        out_file = self.out_dir + self.filename
        if '.npz' not in out_file: out_file += '.npz'
        
        fig_file = self.fig_dir + self.filename
        if '.png' not in fig_file: fig_file += '.png'
        
        
        if (self.B_or_T[0] == 'F') or (self.B_or_T[0] == 'B'):
            B = np.linspace(self.B_r[0], self.B_r[1], self.NB)
            self.out_fileB = out_file.replace('.npz', '_B.npz') if self.B_or_T[0] == 'B' else out_file
            self.fig_fileB = fig_file.replace('.png', '_B.png') if self.B_or_T[0] == 'B' else fig_file
            rho = B_sweep(B, self.TB, self.N_k, self.dt, self.N_tau, self.out_fileB, 
                          self.s_coeff, self.tau_choice, self.tau_coeff, self.in_AFBZ, self.FLEX)
            
        if (self.B_or_T[0] == 'T') or (self.B_or_T[0] == 'B'):
            T = np.linspace(self.T_r[0], self.T_r[1], self.NT)
            self.out_fileT = out_file.replace('.npz', '_T.npz') if self.B_or_T[0] == 'B' else out_file
            self.fig_fileT = fig_file.replace('.png', '_T.png') if self.B_or_T[0] == 'B' else fig_file
            rho = T_sweep(T, self.BT, self.N_k, self.dt, self.N_tau, self.out_fileT,
                          self.s_coeff, self.tau_choice, self.tau_coeff, self.in_AFBZ, self.FLEX)
        
        self.plot_results()