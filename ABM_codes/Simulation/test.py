# test.py - SIMPLIFIED with corrected units for biomass concentration

from fipy import CellVariable, Grid3D, TransientTerm, DiffusionTerm, ImplicitSourceTerm
import numpy as np
import time
import variablevals as params


class FiPyBiofilmSolver:
    def __init__(self):
        self.mesh = Grid3D(nx=params.nx, ny=params.ny, nz=params.nz,
                           dx=params.dx, dy=params.dy, dz=params.dz)
        n = params.nx * params.ny * params.nz

        # Metabolite fields
        self.sucrose = CellVariable(mesh=self.mesh, value=np.zeros(n))
        self.lactate = CellVariable(mesh=self.mesh, value=np.zeros(n))
        self.acetate = CellVariable(mesh=self.mesh, value=np.zeros(n))
        self.propion = CellVariable(mesh=self.mesh, value=np.zeros(n))
        self.eps     = CellVariable(mesh=self.mesh, value=np.zeros(n))

        # Biomass fields (kg/voxel)
        self.sm_mass = CellVariable(mesh=self.mesh, value=np.zeros(n))
        self.vp_mass = CellVariable(mesh=self.mesh, value=np.zeros(n))

        # Boundary conditions
        back = self.mesh.facesBack
        self.sucrose.constrain(params.INITIAL_SUCROSE_TOP, back)

        front = self.mesh.facesFront
        for var in (self.sucrose, self.lactate, self.acetate, self.propion, self.eps):
            var.faceGrad.constrain(0, front)

        from fipy.solvers import DefaultSolver
        self.solver = DefaultSolver(tolerance=params.fipy_solver_tolerance,
                                    iterations=params.fipy_max_iterations)
        self.epsilon = 1e-8

        # Precompute voxel volume in L
        self.voxel_volume_L = (params.dx * 1e-6) * (params.dy * 1e-6) * (params.dz * 1e-6) * 1000.0  # L

    @staticmethod
    def _to_1d(a):
        a = np.asarray(a)
        return np.ravel(a, order='F') if a.ndim == 3 else a

    def step(self, sucrose_field, lactate_field, acetate_field, propionate_field, eps_field,
             sm_mass_init, vp_mass_init, relaxation_steps=None):

        if relaxation_steps is None:
            relaxation_steps = params.relaxationmcs

        # Load values from CC3D
        self.sucrose.value[:] = self._to_1d(sucrose_field)
        self.lactate.value[:] = self._to_1d(lactate_field)
        self.acetate.value[:] = self._to_1d(acetate_field)
        self.propion.value[:] = self._to_1d(propionate_field)
        self.eps.value[:]     = self._to_1d(eps_field)

        self.sm_mass.value[:] = self._to_1d(sm_mass_init).astype(float)  # kg/voxel
        self.vp_mass.value[:] = self._to_1d(vp_mass_init).astype(float)  # kg/voxel

        initial_sm_mass = self.sm_mass.value.copy()
        initial_vp_mass = self.vp_mass.value.copy()

        total_time = params.bacterial_update_interval  # s
        dt_relax   = params.metabolite_timestep        # s

        wall0 = time.time()
        print(f"    FiPy: relaxation_steps={relaxation_steps}, dt={dt_relax:.3f}s")

        for step in range(relaxation_steps):
            t_step = time.time()

            # Convert biomass (kg/voxel) → concentration (kg/L)
            Xsm = self.sm_mass / self.voxel_volume_L
            Xvp = self.vp_mass / self.voxel_volume_L

            # Inhibition factor
            phi = 1.0 / (1.0 + self.lactate/params.K_I_La +
                         self.acetate/params.K_I_Ac + self.propion/params.K_I_Pr)

            # Reaction rates (all in per-second consistent units)
            su_cons = (params.V_max_SM_Su_mass_s * Xsm * phi *
                       self.sucrose / (params.K_M_SM_Su + self.sucrose + self.epsilon))
            la_prod = params.Y_La_Su * su_cons
            la_cons = (params.V_max_VP_La_mass_s * Xvp *
                       self.lactate / (params.K_M_VP_La + self.lactate + self.epsilon))
            ac_prod = (params.Y_Ac_La * params.V_max_VP_Ac_mass_s * Xvp *
                       self.lactate / (params.K_M_VP_Ac + self.lactate + self.epsilon))
            pr_prod = (params.Y_Pr_La * params.V_max_VP_Pr_mass_s * Xvp *
                       self.lactate / (params.K_M_VP_Pr + self.lactate + self.epsilon))

            eps_rate_mass = (params.q_E0_mass_s +
                             params.beta_mass_s * (Xvp / (params.K_beta_mass + Xvp + 1e-18)))
            eps_prod = eps_rate_mass * Xsm

            # PDEs
            eq_su  = (TransientTerm() == DiffusionTerm(params.D_Su)  - su_cons)
            eq_la  = (TransientTerm() == DiffusionTerm(params.D_La)  + la_prod - la_cons)
            eq_ac  = (TransientTerm() == DiffusionTerm(params.D_Ac)  + ac_prod
                      - ImplicitSourceTerm(params.delta_Ac_s))
            eq_pr  = (TransientTerm() == DiffusionTerm(params.D_Pr)  + pr_prod
                      - ImplicitSourceTerm(params.delta_Pr_s))
            eq_eps = (TransientTerm() == DiffusionTerm(params.D_EPS) + eps_prod
                      - ImplicitSourceTerm(params.delta_E_s))

            eq_su.solve(var=self.sucrose, dt=dt_relax, solver=self.solver)
            eq_la.solve(var=self.lactate, dt=dt_relax, solver=self.solver)
            eq_ac.solve(var=self.acetate, dt=dt_relax, solver=self.solver)
            eq_pr.solve(var=self.propion, dt=dt_relax, solver=self.solver)
            eq_eps.solve(var=self.eps,     dt=dt_relax, solver=self.solver)

        # Biomass ODEs
        Xsm = self.sm_mass / self.voxel_volume_L
        Xvp = self.vp_mass / self.voxel_volume_L

        phi = 1.0 / (1.0 + self.lactate/params.K_I_La +
                     self.acetate/params.K_I_Ac + self.propion/params.K_I_Pr)

        sm_growth = (params.mu_max_SM_s * self.sucrose /
                     (params.K_S_SM_Su + self.sucrose + self.epsilon) * phi)
        eps_enh   = params.k_E_s * self.eps
        sm_net    = sm_growth + eps_enh - params.death_SM_s

        vp_growth = (params.mu_max_VP_s * self.lactate /
                     (params.K_S_VP_La + self.lactate + self.epsilon))
        vp_net    = vp_growth - params.death_VP_s

        eq_sm = (TransientTerm() == sm_net * self.sm_mass)
        eq_vp = (TransientTerm() == vp_net * self.vp_mass)
        eq_sm.solve(var=self.sm_mass, dt=total_time, solver=self.solver)
        eq_vp.solve(var=self.vp_mass, dt=total_time, solver=self.solver)

        # Non-negativity enforcement
        for var in (self.sm_mass, self.vp_mass, self.sucrose,
                    self.lactate, self.acetate, self.propion, self.eps):
            var.setValue(np.maximum(var.value, 0))

        sm_mass_change = self.sm_mass.value - initial_sm_mass
        vp_mass_change = self.vp_mass.value - initial_vp_mass

        print(f"    [FiPy] done. total wall {time.time() - wall0:.2f}s")
        return [
            np.array(self.sucrose.value),
            np.array(self.lactate.value),
            np.array(self.acetate.value),
            np.array(self.propion.value),
            np.array(self.eps.value),
            np.array(self.sm_mass.value),
            np.array(self.vp_mass.value),
            np.array(sm_mass_change),
            np.array(vp_mass_change),
        ]


def tester(*args, **kwargs):
    solver = FiPyBiofilmSolver()
    return solver.step(*args, **kwargs)