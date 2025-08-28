# variablevals.py - Parameters + auto time configuration from scientific targets

import numpy as np
import math

# =============================================================================
# 1) USER-INTENT TIME TARGETS (laptop-friendly)
# =============================================================================
TARGET_TOTAL_TIME_MIN = 60.0
BIOLOGY_UPDATE_FRACTION_OF_DOUBLING = 1.0 / 20.0
RELAX_CFL_SAFETY = 0.2
BASE_S_MCS = 1.0
RELAX_STEPS_MIN = 5
RELAX_STEPS_MAX = 10 # lowered from 20 to speed FiPy
DOUBLING_TIME_MIN_SM = None
DOUBLING_TIME_MIN_VP = None

# =============================================================================
# 2) GRID & DOMAIN PARAMETERS
# =============================================================================
S_MUTANS = 1
V_PARVULA = 2

nx = 50
ny = 50
nz = 14
dx = 1.0  # µm
dy = 1.0  # µm
dz = 1.0  # µm

# =============================================================================
# 3) BIOCHEMICAL KINETIC PARAMETERS (per hour unless noted)
# =============================================================================
# Metabolite concentration unit: mM (mmol/L)
# EPS concentration unit: mg/L

# --- S. mutans ---
mu_max_SM = 5.0       # 1/hr
K_S_SM_Su = 1.0       # mM
V_max_SM_Su = 0.1    # mmol/(kg·hr)
K_M_SM_Su = 0.12      # mM
Y_La_Su = 0.3         # mol/mol
death_SM = 0.0001      # 1/hr

# --- V. parvula ---
mu_max_VP = 4.0       # 1/hr
K_S_VP_La = 2.42      # mM
V_max_VP_La = 0.1    # mmol/(kg·hr)
K_M_VP_La = 0.1       # mM
death_VP = 0.00005     # 1/hr

# --- Acid production / inhibition ---
V_max_VP_Ac = 0.01    # mmol/(kg·hr)
K_M_VP_Ac   = 0.5     # mM
delta_Ac    = 0.005   # 1/hr
Y_Ac_La     = 0.5     # mol/mol

V_max_VP_Pr = 0.01    # mmol/(kg·hr)
K_M_VP_Pr   = 0.5     # mM
delta_Pr    = 0.005   # 1/hr
Y_Pr_La     = 0.3     # mol/mol

K_I_La = 25.0         # mM
K_I_Ac = 20.0         # mM
K_I_Pr = 20.0         # mM

# --- EPS ---
q_E0   = 0.01         # mg/(kg·hr)
beta   = 0.5          # mg/(kg·hr)
K_beta_mass = 1e-18   # kg/µm³ (tune)
delta_E = 0.001       # 1/hr
k_E    = 0.01         # 1/(hr·(mg/L))

# =============================================================================
# 4) DIFFUSION COEFFICIENTS (µm²/s)
# =============================================================================
D_Su  = 0.1
D_La  = 0.05
D_Ac  = 0.05
D_Pr  = 0.05
D_EPS = 0.01
D_MAX = max(D_Su, D_La, D_Ac, D_Pr, D_EPS)

# -----------------------------------------------------------------------------
# 5) MASS–VOLUME COUPLING PARAMETERS
# -----------------------------------------------------------------------------
VOXEL_LENGTH_UM = 1.0
DX = VOXEL_LENGTH_UM
DY = VOXEL_LENGTH_UM
DZ = VOXEL_LENGTH_UM

MASS_PER_BACTERIUM = 1.1e-15  # kg

VOXEL_VOLUME_M3 = (VOXEL_LENGTH_UM * 1e-6)**3
VOXEL_VOLUME_UM3 = VOXEL_LENGTH_UM**3

BACTERIAL_DENSITY = MASS_PER_BACTERIUM / VOXEL_VOLUME_M3

SINGLE_BACTERIUM_VOLUME_M3 = MASS_PER_BACTERIUM / BACTERIAL_DENSITY
SINGLE_BACTERIUM_VOLUME_UM3 = SINGLE_BACTERIUM_VOLUME_M3 * 1e18

LUMP_FACTOR = 1.0
LUMP_AFFECTS_GEOMETRY = False

DENSITY_GEOM = BACTERIAL_DENSITY * (LUMP_FACTOR if LUMP_AFFECTS_GEOMETRY else 1.0)

MIN_CELL_VOLUME = 2
MAX_CELL_VOLUME = 64
DIVISION_VOLUME = 20  # unused (reference)

# Mass-based division - REASONABLE THRESHOLDS
DIVISION_MASS_IN_BACTERIA = 10.0   # Divide at 10 bacteria worth
DAUGHTER_MASS_IN_BACTERIA = 5.0     # Daughters get 5 bacteria worth
DIVISION_MASS_THRESHOLD_KG = DIVISION_MASS_IN_BACTERIA * MASS_PER_BACTERIUM
DAUGHTER_CELL_MASS_KG = DAUGHTER_MASS_IN_BACTERIA * MASS_PER_BACTERIUM

MASS_PER_VOXEL = DENSITY_GEOM * VOXEL_VOLUME_M3

def mass_to_volume(mass_kg):
    if mass_kg <= 0:
        return MIN_CELL_VOLUME
    volume_m3 = mass_kg / DENSITY_GEOM
    volume_um3 = volume_m3 * 1e18
    volume_vox = volume_um3 / VOXEL_VOLUME_UM3
    return max(MIN_CELL_VOLUME, min(MAX_CELL_VOLUME, volume_vox))

def volume_to_mass(volume_voxels):
    volume_um3 = volume_voxels * VOXEL_VOLUME_UM3
    volume_m3  = volume_um3 * 1e-18
    return volume_m3 * DENSITY_GEOM

# =============================================================================
# 5b) CELL SEEDING PARAMETERS
# =============================================================================
INITIAL_SEED_SIZE   = 2
INITIAL_CELL_VOLUME = 4.0

# =============================================================================
# 6) INITIAL CONDITIONS (mM for metabolites, mg/L for EPS)
# =============================================================================
INITIAL_SM_MIN = 8
INITIAL_SM_MAX = 12
INITIAL_VP_MIN = 8
INITIAL_VP_MAX = 12

EXPERIMENT_TYPE = "separated"  # Options: "mixed", "separated", "sm_only", "vp_only"

Y_Ac_La = 0.5
Y_Pr_La = 0.3

ACETATE_GROWTH_ENHANCEMENT = 0.2
PROPIONATE_GROWTH_ENHANCEMENT = 0.1

EPS_ENHANCEMENT_BY_ACETATE = 1.5
EPS_ENHANCEMENT_BY_PROPIONATE = 1.2

INITIAL_SEPARATION_DISTANCE = 10
CLUSTER_RADIUS = 5

INITIAL_SUCROSE_BOTTOM = 0.0
INITIAL_SUCROSE_TOP    = 5.0

initial_lactate    = 0.0
initial_acetate    = 0.0
initial_propionate = 0.0
initial_eps        = 0.0

# =============================================================================
# 7) AUTO TIME CONFIGURATION
# =============================================================================
TARGET_TOTAL_TIME_S = float(TARGET_TOTAL_TIME_MIN) * 60.0

def doubling_time_sec_from_mu(mu_per_hr):
    if mu_per_hr <= 0:
        return math.inf
    return (math.log(2.0) / mu_per_hr) * 3600.0

DT_SM_S = (DOUBLING_TIME_MIN_SM * 60.0) if DOUBLING_TIME_MIN_SM else doubling_time_sec_from_mu(mu_max_SM)
DT_VP_S = (DOUBLING_TIME_MIN_VP * 60.0) if DOUBLING_TIME_MIN_VP else doubling_time_sec_from_mu(mu_max_VP)

FASTEST_DT_S = min(DT_SM_S, DT_VP_S)
if not math.isfinite(FASTEST_DT_S) or FASTEST_DT_S <= 0:
    FASTEST_DT_S = 3600.0

delta_t_bio_target = FASTEST_DT_S * BIOLOGY_UPDATE_FRACTION_OF_DOUBLING

CFL_LIMIT_S = (dx ** 2) / (2.0 * max(D_MAX, 1e-12))
dt_relax_target = max(1e-6, min(CFL_LIMIT_S * RELAX_CFL_SAFETY, delta_t_bio_target / 20.0))

s_mcs = float(BASE_S_MCS)
steppable_frequency = max(1, int(round(delta_t_bio_target / s_mcs)))
bacterial_update_interval = steppable_frequency * s_mcs

relaxationmcs = int(np.clip(int(math.ceil(bacterial_update_interval / dt_relax_target)),
                            RELAX_STEPS_MIN, RELAX_STEPS_MAX))
metabolite_timestep = bacterial_update_interval / relaxationmcs
RECOMMENDED_STEPS = int(math.ceil(TARGET_TOTAL_TIME_S / s_mcs))
RELAXATION_RATIO = steppable_frequency / max(1, relaxationmcs)

h_mcs = s_mcs / 3600.0  # hours per MCS

# =============================================================================
# 8) UNITS & PER-SECOND CONVERSIONS
# =============================================================================
H2S = 1.0 / 3600.0
UM3_PER_L = 1e15

mu_max_SM_s = mu_max_SM * H2S
mu_max_VP_s = mu_max_VP * H2S
death_SM_s  = death_SM  * H2S
death_VP_s  = death_VP  * H2S

V_max_SM_Su_mass_s = V_max_SM_Su * H2S
V_max_VP_La_mass_s = V_max_VP_La * H2S
V_max_VP_Ac_mass_s = V_max_VP_Ac * H2S
V_max_VP_Pr_mass_s = V_max_VP_Pr * H2S

delta_Ac_s = delta_Ac * H2S
delta_Pr_s = delta_Pr * H2S
delta_E_s  = delta_E  * H2S

q_E0_mass_s = q_E0 * H2S
beta_mass_s = beta * H2S
k_E_s = k_E * H2S

# =============================================================================
# 9) FIPY SOLVER CONFIG
# =============================================================================
fipy_solver_tolerance = 1e-5
fipy_max_iterations   = 50   # lowered from 1000 to speed up
relaxation_convergence_criterion = 0.01

# =============================================================================
# 10) PERFORMANCE + VERBOSITY
# =============================================================================
VERBOSE_MCS_LOG = True       # print every MCS heartbeat
VERBOSE_RELAXATION = True    # print FiPy progress
VERBOSE_SUBSTEPS = False      # print each substep inside a relaxation step
RELAXATION_LOG_EVERY = 10     # print every relaxation step
PLOT_FREQUENCY = 50
BIOLOGY_START_MCS = 0
FIPY_FREQUENCY = 100       # FiPy runs every 500 MCS
HEARTBEAT_EVERY_MCS = 1      # heartbeat every MCS

# =============================================================================
# 11) PRINT CONFIG
# =============================================================================
print("=" * 60)
print("BIOFILM SIMULATION CONFIGURATION (auto time setup)")
print("=" * 60)
print(f"Domain: {nx}x{ny}x{nz} grid ({nx*dx}x{ny*dy}x{nz*dz} µm)")
print(f"D_max: {D_MAX:.3f} µm²/s; CFL limit: {CFL_LIMIT_S:.3f} s")
print("-" * 60)
print(f"Target total story time: {TARGET_TOTAL_TIME_MIN:.1f} min ({TARGET_TOTAL_TIME_S:.1f} s)")
print(f"Base seconds per MCS: s_mcs = {s_mcs:.4f} s")
print(f"Biology update interval target: {delta_t_bio_target:.3f} s "
      f"(fraction {BIOLOGY_UPDATE_FRACTION_OF_DOUBLING:.4f} of fastest doubling)")
print(f"Chosen steppable_frequency: {steppable_frequency} MCS "
      f"(actual Δt_bio = {bacterial_update_interval:.4f} s)")
print(f"Relaxation inner steps: {relaxationmcs} (dt_relax = {metabolite_timestep:.4f} s)")
print(f"Relaxation/Chemistry safety fraction of CFL: {RELAX_CFL_SAFETY:.2f}")
print("-" * 60)
print(f"Recommended XML Steps: {RECOMMENDED_STEPS} (set <Steps>{RECOMMENDED_STEPS}</Steps> in XML)")
print(f"Timescale separation: steppable_frequency / relaxationmcs = {RELAXATION_RATIO:.2f}:1")
print("-" * 60)
print(f"Single bacterium volume: {SINGLE_BACTERIUM_VOLUME_UM3:.1f} µm³")
print(f"Initial cell volume: {INITIAL_CELL_VOLUME} voxels")
print(f"Initial cell counts: SM=[{INITIAL_SM_MIN},{INITIAL_SM_MAX}], VP=[{INITIAL_VP_MIN},{INITIAL_VP_MAX}]")
print(f"Initial lactate (mM): {initial_lactate}")
print("=" * 60)