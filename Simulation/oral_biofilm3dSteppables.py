# oral_biofilm3dSteppables.py - Mass-driven volume update + pixel-based mass mapping

import os
import numpy as np
import random
import time
import math
from cc3d.core.PySteppables import *
import variablevals as params
from test import FiPyBiofilmSolver

class oral_biofilm3dSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.stats = {
            'field_extraction': [], 'mass_extraction': [], 'fipy_solver': [],
            'volume_update': [], 'field_update': [], 'plotting': [], 'total_step': []
        }
        self.bio_frequency = params.steppable_frequency
        self.solver = None
        self.fipy_frequency = int(getattr(params, 'FIPY_FREQUENCY', 1000))
        self.heartbeat_every = int(getattr(params, 'HEARTBEAT_EVERY_MCS', 1))
        self.log_path = None

    def _log(self, msg):
        # Print and append to run.log (flush immediately)
        print(msg, flush=True)
        if self.log_path:
            with open(self.log_path, "a", encoding="utf-8") as f:
                f.write(str(msg) + "\n")
                f.flush()

    def start(self):
        random.seed(42); np.random.seed(42)
        self.log_path = os.path.join(self.output_dir, "run.log")
        self._log("=== STARTING BIOFILM SIMULATION (mass-driven volumes) ===")

        t0 = time.time()
        self._initialize_fields()
        self._setup_plots()
        self.create_bacterial_cells()
        self._initialize_cell_masses()

        self.solver = FiPyBiofilmSolver()

        self._log(f"Initialization completed in {time.time() - t0:.2f} s")
        self._diagnose_sucrose_field(0)

    def create_bacterial_cells(self):
        t0 = time.time()
        n_smutans = random.randint(params.INITIAL_SM_MIN, params.INITIAL_SM_MAX)
        n_vparvula = random.randint(params.INITIAL_VP_MIN, params.INITIAL_VP_MAX)

        experiment_type = getattr(params, 'EXPERIMENT_TYPE', 'mixed')
        self._log(f"=== EXPERIMENT TYPE: {experiment_type} ===")
        self._log(f"Seeding {n_smutans} S. mutans, {n_vparvula} V. parvula")

        seed_len = int(max(1, params.INITIAL_SEED_SIZE))
        z_level = 0

        placed_sm = 0
        placed_vp = 0

        if experiment_type == "sm_only":
            n_vparvula = 0
            placed_sm = self._place_cells_in_circle(params.S_MUTANS, n_smutans,
                                                    self.dim.x//2, self.dim.y//2,
                                                    min(self.dim.x, self.dim.y)//4, z_level)
        elif experiment_type == "vp_only":
            n_smutans = 0
            placed_vp = self._place_cells_in_circle(params.V_PARVULA, n_vparvula,
                                                    self.dim.x//2, self.dim.y//2,
                                                    min(self.dim.x, self.dim.y)//4, z_level)
        elif experiment_type == "separated":
            separation = getattr(params, 'INITIAL_SEPARATION_DISTANCE', 10)
            cluster_r = getattr(params, 'CLUSTER_RADIUS', 5)
            placed_sm = self._place_cells_in_circle(params.S_MUTANS, n_smutans,
                                                    self.dim.x//2 - separation//2, self.dim.y//2,
                                                    cluster_r, z_level)
            placed_vp = self._place_cells_in_circle(params.V_PARVULA, n_vparvula,
                                                    self.dim.x//2 + separation//2, self.dim.y//2,
                                                    cluster_r, z_level)
        elif experiment_type == "mixed":
            placed_sm, placed_vp = self._place_cells_mixed(n_smutans, n_vparvula, z_level)
        elif experiment_type == "layered":
            placed_sm = self._place_cells_in_circle(params.S_MUTANS, n_smutans,
                                                    self.dim.x//2, self.dim.y//2,
                                                    min(self.dim.x, self.dim.y)//4, 0)
            placed_vp = self._place_cells_in_circle(params.V_PARVULA, n_vparvula,
                                                    self.dim.x//2, self.dim.y//2,
                                                    min(self.dim.x, self.dim.y)//4, 2)

        self._log(f"Cell seeding completed in {time.time() - t0:.2f} s")
        self._log(f"Placed: SM={placed_sm}, VP={placed_vp}, Total={placed_sm+placed_vp}")
        self._log(f"Configuration: {experiment_type}")

    def _place_cells_in_circle(self, cell_type, num_cells, center_x, center_y, radius, z_level):
        placed = 0
        seed_len = int(max(1, params.INITIAL_SEED_SIZE))
        occupied = []
        for _ in range(num_cells):
            for _attempt in range(100):
                angle = random.uniform(0, 2 * np.pi)
                r = random.uniform(0, radius)
                x = int(center_x + r * math.cos(angle))
                y = int(center_y + r * math.sin(angle))
                z = z_level

                if not (seed_len <= x < self.dim.x - seed_len and seed_len <= y < self.dim.y - seed_len):
                    continue

                too_close = False
                for (ox, oy, oz) in occupied:
                    if ((x - ox)**2 + (y - oy)**2 + (z - oz)**2)**0.5 < seed_len * 1.5:
                        too_close = True
                        break

                if not too_close:
                    cell = self.new_cell(cell_type)
                    cell.targetVolume = float(params.INITIAL_CELL_VOLUME)
                    cell.lambdaVolume = 2.0

                    initial_mass = params.volume_to_mass(cell.targetVolume)
                    cell.dict['mass'] = initial_mass
                    cell.dict['mass_at_last_volume_update'] = initial_mass
                    cell.dict['mass_change'] = 0.0
                    cell.dict['volume_change'] = 0.0
                    cell.dict['age'] = 0
                    cell.dict['generation'] = 0

                    self.cell_field[x:x+seed_len, y:y+seed_len, z:z+seed_len] = cell
                    occupied.append((x, y, z))
                    placed += 1
                    break
        return placed

    def _place_cells_mixed(self, n_sm, n_vp, z_level):
        placed_sm = 0
        placed_vp = 0
        seed_len = int(max(1, params.INITIAL_SEED_SIZE))
        center_x = self.dim.x // 2
        center_y = self.dim.y // 2
        radius = min(self.dim.x, self.dim.y) // 4

        total = n_sm + n_vp
        cells_to_place = []
        for i in range(total):
            if i % 2 == 0 and placed_sm < n_sm:
                cells_to_place.append(params.S_MUTANS); placed_sm += 1
            elif placed_vp < n_vp:
                cells_to_place.append(params.V_PARVULA); placed_vp += 1
            elif placed_sm < n_sm:
                cells_to_place.append(params.S_MUTANS); placed_sm += 1

        placed_sm = 0
        placed_vp = 0
        occupied = []
        for cell_type in cells_to_place:
            for _attempt in range(100):
                angle = random.uniform(0, 2 * np.pi)
                r = random.uniform(0, radius)
                x = int(center_x + r * math.cos(angle))
                y = int(center_y + r * math.sin(angle))
                z = z_level
                if not (seed_len <= x < self.dim.x - seed_len and seed_len <= y < self.dim.y - seed_len):
                    continue
                too_close = False
                for (ox, oy, oz) in occupied:
                    if ((x - ox)**2 + (y - oy)**2 + (z - oz)**2)**0.5 < seed_len * 1.5:
                        too_close = True
                        break
                if not too_close:
                    cell = self.new_cell(cell_type)
                    cell.targetVolume = float(params.INITIAL_CELL_VOLUME)
                    cell.lambdaVolume = 2.0
                    initial_mass = params.volume_to_mass(cell.targetVolume)
                    cell.dict['mass'] = initial_mass
                    cell.dict['mass_at_last_volume_update'] = initial_mass
                    cell.dict['mass_change'] = 0.0
                    cell.dict['volume_change'] = 0.0
                    cell.dict['age'] = 0
                    cell.dict['generation'] = 0
                    self.cell_field[x:x+seed_len, y:y+seed_len, z:z+seed_len] = cell
                    occupied.append((x, y, z))
                    if cell_type == params.S_MUTANS: placed_sm += 1
                    else: placed_vp += 1
                    break
        return placed_sm, placed_vp

    def _initialize_cell_masses(self):
        for cell in self.cell_list:
            if cell.type in (params.S_MUTANS, params.V_PARVULA):
                if 'mass' not in cell.dict:
                    cell.dict['mass'] = params.volume_to_mass(cell.targetVolume)
                cell.dict.setdefault('mass_at_last_volume_update', cell.dict['mass'])
                cell.dict.setdefault('mass_change', 0.0)
                cell.dict.setdefault('volume_change', 0.0)
                cell.dict.setdefault('age', 0)
                cell.dict.setdefault('generation', 0)

    def step(self, mcs):
        step_t0 = time.time()

        # Heartbeat every MCS
        if self.heartbeat_every > 0 and (mcs % self.heartbeat_every == 0):
            self._log(f"[MCS {mcs}] alive; cells={len(self.cell_list)}")

        # Plots
        plot_freq = getattr(params, 'PLOT_FREQUENCY', 50)
        if mcs == 0 or mcs % plot_freq == 0:
            t = time.time()
            self._update_plots(mcs)
            self.stats['plotting'].append(time.time() - t)

        # Occasional detailed dump
        if mcs % 500 == 0 and mcs > 0:
            self._debug_cells(mcs)

        bio_start = getattr(params, 'BIOLOGY_START_MCS', 0)
        do_bio = (mcs % self.bio_frequency == 0) and (mcs >= bio_start)

        if not do_bio:
            if getattr(params, 'VERBOSE_MCS_LOG', False):
                self._log(f"[MCS {mcs}] skipping biology (start={bio_start}, freq={self.bio_frequency})")
            self.stats['total_step'].append(time.time() - step_t0)
            return

        # Biology update phases
        self._log(f"\n=== BIOLOGY UPDATE at MCS {mcs} ===")
        try:
            self._log("  1. Extracting field data...")
            t = time.time()
            current_fields = self._extract_field_data()
            if current_fields is None:
                self._log("  ERROR: field extraction failed")
                return
            self._log(f"     Sucrose avg: {np.mean(current_fields[0]):.6f}")
            self._log(f"     Lactate avg: {np.mean(current_fields[1]):.6f}")
            self.stats['field_extraction'].append(time.time() - t)

            self._log("  2. Extracting cell masses...")
            t = time.time()
            sm_mass_field, vp_mass_field = self._extract_cell_mass_fields()
            self._log(f"     SM total mass: {np.sum(sm_mass_field):.3e}")
            self._log(f"     VP total mass: {np.sum(vp_mass_field):.3e}")
            self.stats['mass_extraction'].append(time.time() - t)

            if mcs % self.fipy_frequency == 0 and mcs > 0:
                self._log("  3. Running FiPy solver...")
                t = time.time()
                if self.solver is None:
                    self._log("     Solver missing, creating...")
                    self.solver = FiPyBiofilmSolver()

                fipy_t0 = time.time()
                updated_fields = self._run_relaxation_solver_with_mass(
                    current_fields, sm_mass_field, vp_mass_field
                )
                self._log(f"     FiPy wall: {time.time() - fipy_t0:.2f}s")

                if not isinstance(updated_fields, (list, tuple)) or len(updated_fields) < 9:
                    self._log("     ERROR: FiPy returned incomplete data; skipping updates.")
                else:
                    self._log(f"     Returned {len(updated_fields)} fields")
                    self._log(f"     Σ dM_SM={np.sum(updated_fields[7]):.3e}, Σ dM_VP={np.sum(updated_fields[8]):.3e}")
                    self.stats['fipy_solver'].append(time.time() - t)

                    self._log("  4. Updating CC3D fields...")
                    t = time.time()
                    self._update_fields(updated_fields)
                    self.stats['field_update'].append(time.time() - t)

                    self._log("  5. Updating cell volumes...")
                    t = time.time()
                    self._update_cell_volumes_from_mass_changes(updated_fields)
                    self.stats['volume_update'].append(time.time() - t)
            else:
                self._log(f"  3. Skipping FiPy (runs every {self.fipy_frequency} MCS)")
                # Do nothing - let FiPy handle all growth naturally when it runs

            self._print_biofilm_status_with_volumes(mcs)
            self.stats['total_step'].append(time.time() - step_t0)

            if len(self.stats['total_step']) % 10 == 0:
                self._print_timing_stats()

        except Exception as e:
            self._log(f"ERROR in step {mcs}: {e}")
            import traceback
            traceback.print_exc()

    def _extract_cell_mass_fields(self):
        sm = np.zeros((params.nx, params.ny, params.nz), dtype=float)
        vp = np.zeros_like(sm)
        for cell in self.cell_list:
            if cell.type not in (params.S_MUTANS, params.V_PARVULA):
                continue
            mass = cell.dict.get('mass', params.volume_to_mass(cell.targetVolume))
            plist = list(self.get_cell_pixel_list(cell))
            if not plist:
                continue
            share = mass / float(len(plist))
            target = sm if cell.type == params.S_MUTANS else vp
            for ptd in plist:
                pt = ptd.pixel
                x, y, z = pt.x, pt.y, pt.z
                if 0 <= x < params.nx and 0 <= y < params.ny and 0 <= z < params.nz:
                    target[x, y, z] += share
        return sm, vp

    def _update_cell_volumes_from_mass_changes(self, updated_fields):
        if not updated_fields or len(updated_fields) < 9:
            self._log("Warning: No mass deltas from solver.")
            return

        sm_d = np.array(updated_fields[7]).reshape((params.nx, params.ny, params.nz), order='F')
        vp_d = np.array(updated_fields[8]).reshape((params.nx, params.ny, params.nz), order='F')

        grew = shrank = divided = died = 0
        MASS_PER_10_BACTERIA = 10 * params.MASS_PER_BACTERIUM

        for cell in list(self.cell_list):
            if not cell:
                continue
            if cell.type not in (params.S_MUTANS, params.V_PARVULA):
                continue

            old_mass = cell.dict.get('mass', 0.0)
            plist = list(self.get_cell_pixel_list(cell))
            dm = 0.0
            if plist:
                d_field = sm_d if cell.type == params.S_MUTANS else vp_d
                for ptd in plist:
                    pt = ptd.pixel
                    if 0 <= pt.x < params.nx and 0 <= pt.y < params.ny and 0 <= pt.z < params.nz:
                        dm += d_field[pt.x, pt.y, pt.z]

            new_mass = max(0.0, old_mass + dm)
            cell.dict['mass'] = new_mass
            cell.dict['mass_change'] = dm

            mass_at_last_update = cell.dict.get('mass_at_last_volume_update', old_mass)
            mass_accumulated = new_mass - mass_at_last_update
            volume_change = int(mass_accumulated / MASS_PER_10_BACTERIA)

            if volume_change != 0:
                old_vol = float(cell.targetVolume)
                new_vol = old_vol + volume_change
                new_vol = max(params.MIN_CELL_VOLUME, min(params.MAX_CELL_VOLUME, new_vol))
                cell.targetVolume = new_vol
                cell.dict['volume_change'] = new_vol - old_vol
                mass_used = volume_change * MASS_PER_10_BACTERIA
                cell.dict['mass_at_last_volume_update'] = mass_at_last_update + mass_used
                grew += 1 if new_vol > old_vol else 0
                shrank += 1 if new_vol < old_vol else 0

            cell.dict['age'] = cell.dict.get('age', 0) + 1

            if cell.dict['mass'] >= params.DIVISION_MASS_THRESHOLD_KG:
                if self._divide_cell(cell):
                    divided += 1
            elif cell.targetVolume < params.MIN_CELL_VOLUME:
                self.delete_cell(cell)
                died += 1

        self._log(f"  Cells grew: {grew}, shrank: {shrank}, divided: {divided}, died: {died}")

    def _divide_cell(self, parent_cell):
        daughter = self.new_cell(parent_cell.type)
        if not self.move_cell_into_new_space(daughter, parent_cell):
            self._log(f"  Division failed for cell {parent_cell.id}, no space found.")
            self.delete_cell(daughter)
            return False

        target_mass = params.DAUGHTER_CELL_MASS_KG
        target_volume = params.mass_to_volume(target_mass)

        self._log(f"  Cell {parent_cell.id} (mass {parent_cell.dict['mass']:.2e}) is dividing.")

        parent_cell.targetVolume = target_volume
        parent_cell.lambdaVolume = 2.0
        parent_cell.dict['mass'] = target_mass
        parent_cell.dict['mass_at_last_volume_update'] = target_mass
        parent_cell.dict['age'] = 0
        parent_cell.dict['generation'] = parent_cell.dict.get('generation', 0) + 1

        daughter.targetVolume = target_volume
        daughter.lambdaVolume = 2.0
        daughter.dict['mass'] = target_mass
        daughter.dict['mass_at_last_volume_update'] = target_mass
        daughter.dict['age'] = 0
        daughter.dict['generation'] = parent_cell.dict.get('generation', 0)
        daughter.dict['mass_change'] = 0.0
        daughter.dict['volume_change'] = 0.0

        self._log(f"  -> Parent {parent_cell.id} and Daughter {daughter.id} reset to mass {target_mass:.2e} kg.")
        return True

    def move_cell_into_new_space(self, daughter, parent_cell):
        plist = list(self.get_cell_pixel_list(parent_cell))
        if not plist:
            return False
        seed_pt = random.choice(plist).pixel
        n_parent = len(plist)
        n_move_max = max(1, min(8, n_parent // 4))
        moved = 0
        for dz in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if moved >= n_move_max:
                        break
                    x = seed_pt.x + dx
                    y = seed_pt.y + dy
                    z = seed_pt.z + dz
                    if not (0 <= x < self.dim.x and 0 <= y < self.dim.y and 0 <= z < self.dim.z):
                        continue
                    current = self.cell_field[x, y, z]
                    if current is parent_cell:
                        self.cell_field[x, y, z] = daughter
                        moved += 1
                if moved >= n_move_max:
                    break
            if moved >= n_move_max:
                break
        if moved == 0:
            return False
        if len(list(self.get_cell_pixel_list(parent_cell))) == 0:
            # revert if we accidentally emptied the parent
            for dz in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dx in (-1, 0, 1):
                        x = seed_pt.x + dx
                        y = seed_pt.y + dy
                        z = seed_pt.z + dz
                        if 0 <= x < self.dim.x and 0 <= y < self.dim.y and 0 <= z < self.dim.z:
                            if self.cell_field[x, y, z] is daughter:
                                self.cell_field[x, y, z] = parent_cell
            return False
        return True

    def _run_relaxation_solver_with_mass(self, current_fields, sm_mass_field, vp_mass_field):
        self._log(f"  Running FiPy solver with {params.relaxationmcs} inner steps...")
        updated = self.solver.step(
            current_fields[0], current_fields[1],
            current_fields[2], current_fields[3], current_fields[4],
            sm_mass_field, vp_mass_field,
            relaxation_steps=params.relaxationmcs
        )
        s_before = current_fields[0]
        s_after = np.array(updated[0]).reshape((params.nx, params.ny, params.nz), order='F')
        self._verify_relaxation_convergence(s_before, s_after, self.mcs)
        return updated

    def _print_timing_stats(self):
        self._log("\n=== TIMING STATS (avg of last 10, s) ===")
        for k, vals in self.stats.items():
            if vals:
                self._log(f"  {k:20s}: {np.mean(vals[-10:]):.3f}")
        self._log("=" * 40)

    def _count_cells(self):
        sm = vp = 0
        for cell in self.cell_list:
            if cell.type == params.S_MUTANS: sm += 1
            elif cell.type == params.V_PARVULA: vp += 1
        return {'S_mutans': sm, 'V_parvula': vp, 'total': sm+vp}

    def _initialize_fields(self):
        self._log("Initializing fields with sucrose gradient...")
        for x in range(params.nx):
            for y in range(params.ny):
                for z in range(params.nz):
                    grad = params.INITIAL_SUCROSE_BOTTOM + (
                        (params.INITIAL_SUCROSE_TOP - params.INITIAL_SUCROSE_BOTTOM) *
                        z / max(1, params.nz - 1)
                    )
                    self.field.sucrose[x, y, z] = float(grad)
                    self.field.lactate[x, y, z] = params.initial_lactate
                    self.field.acetate[x, y, z] = params.initial_acetate
                    self.field.propionate[x, y, z] = params.initial_propionate
                    self.field.eps[x, y, z] = params.initial_eps

    def _extract_field_data(self):
        try:
            su = np.zeros((params.nx, params.ny, params.nz))
            la = np.zeros_like(su)
            ac = np.zeros_like(su)
            pr = np.zeros_like(su)
            ep = np.zeros_like(su)
            for x in range(params.nx):
                for y in range(params.ny):
                    for z in range(params.nz):
                        su[x, y, z] = float(self.field.sucrose[x, y, z])
                        la[x, y, z] = float(self.field.lactate[x, y, z])
                        ac[x, y, z] = float(self.field.acetate[x, y, z])
                        pr[x, y, z] = float(self.field.propionate[x, y, z])
                        ep[x, y, z] = float(self.field.eps[x, y, z])
            return [su, la, ac, pr, ep]
        except Exception as e:
            self._log(f"Error extracting fields: {e}")
            return None

    def _update_fields(self, fields):
        nx, ny, nz = self.dim.x, self.dim.y, self.dim.z
        try:
            su = np.array(fields[0]).reshape((nx, ny, nz), order='F')
            la = np.array(fields[1]).reshape((nx, ny, nz), order='F')
            ac = np.array(fields[2]).reshape((nx, ny, nz), order='F')
            pr = np.array(fields[3]).reshape((nx, ny, nz), order='F')
            ep = np.array(fields[4]).reshape((nx, ny, nz), order='F')
            
            # Check for overflow/invalid values and cap them
            if np.any(np.isnan(su)) or np.any(np.isinf(su)) or np.max(np.abs(su)) > 1e10:
                self._log(f"WARNING: Invalid sucrose values detected, capping")
                su = np.clip(su, 0, 100)
            if np.any(np.isnan(la)) or np.any(np.isinf(la)) or np.max(np.abs(la)) > 1e10:
                self._log(f"WARNING: Invalid lactate values detected, capping")
                la = np.clip(la, 0, 100)
                
            for z in range(nz):
                for y in range(ny):
                    for x in range(nx):
                        # Cap all values to reasonable ranges
                        self.field.sucrose[x, y, z] = float(np.clip(su[x, y, z], 0.0, 100.0))
                        self.field.lactate[x, y, z] = float(np.clip(la[x, y, z], 0.0, 100.0))
                        self.field.acetate[x, y, z] = float(np.clip(ac[x, y, z], 0.0, 100.0))
                        self.field.propionate[x, y, z] = float(np.clip(pr[x, y, z], 0.0, 100.0))
                        self.field.eps[x, y, z] = float(np.clip(ep[x, y, z], 0.0, 1000.0))
        except Exception as e:
            self._log(f"Error in _update_fields: {e}")

    def _verify_relaxation_convergence(self, before, after, mcs):
        rel = np.abs(after - before) / (before + 1e-8)
        mx = np.max(rel)
        if mx < params.relaxation_convergence_criterion:
            self._log(f"  Relaxation converged (max rel change {mx:.4f})")
            return True
        return False

    def _get_field_array(self, field):
        arr = np.zeros((self.dim.x, self.dim.y, self.dim.z))
        for z in range(self.dim.z):
            for y in range(self.dim.y):
                for x in range(self.dim.x):
                    arr[x, y, z] = field[x, y, z]
        return arr

    def _debug_cells(self, mcs):
        counts = self._count_cells()
        self._log(f"\n=== CELL DEBUG MCS {mcs} ===")
        self._log(f"S. mutans: {counts['S_mutans']}, V. parvula: {counts['V_parvula']}, Total: {counts['total']}")
        sm_vol = []
        sm_mass = []
        vp_vol = []
        vp_mass = []
        for c in self.cell_list:
            if c.type == params.S_MUTANS:
                sm_vol.append(c.targetVolume)
                sm_mass.append(c.dict.get('mass', 0.0))
            elif c.type == params.V_PARVULA:
                vp_vol.append(c.targetVolume)
                vp_mass.append(c.dict.get('mass', 0.0))
        if sm_vol:
            self._log(f"  SM volume min/avg/max: {min(sm_vol):.1f}/{np.mean(sm_vol):.1f}/{max(sm_vol):.1f}")
            self._log(f"  SM mass   min/avg/max: {min(sm_mass):.2e}/{np.mean(sm_mass):.2e}/{max(sm_mass):.2e}")
        if vp_vol:
            self._log(f"  VP volume min/avg/max: {min(vp_vol):.1f}/{np.mean(vp_vol):.1f}/{max(vp_vol):.1f}")
            self._log(f"  VP mass   min/avg/max: {min(vp_mass):.2e}/{np.mean(vp_mass):.2e}/{max(vp_mass):.2e}")

    def _diagnose_sucrose_field(self, mcs):
        self._log(f"\n=== SUCROSE FIELD DIAG at MCS {mcs} ===")
        for z in [0, 7, 13]:
            vals = []
            for x in range(0, params.nx, 10):
                for y in range(0, params.ny, 10):
                    vals.append(float(self.field.sucrose[x, y, z]))
            if vals:
                self._log(f"  z={z:2d}: avg={np.mean(vals):.6f}")

    def _print_biofilm_status_with_volumes(self, mcs):
        counts = self._count_cells()
        sm_cells = [c for c in self.cell_list if c.type == params.S_MUTANS]
        vp_cells = [c for c in self.cell_list if c.type == params.V_PARVULA]

        def stats(cells):
            if not cells:
                return dict(count=0, avg_volume=0, total_volume=0,
                            avg_mass=0, total_mass=0, avg_age=0, max_generation=0)
            vols = [c.targetVolume for c in cells]
            masses = [c.dict.get('mass', 0.0) for c in cells]
            ages = [c.dict.get('age', 0) for c in cells]
            gens = [c.dict.get('generation', 0) for c in cells]
            return dict(
                count=len(cells),
                avg_volume=np.mean(vols), total_volume=np.sum(vols),
                avg_mass=np.mean(masses), total_mass=np.sum(masses),
                avg_age=np.mean(ages), max_generation=max(gens) if gens else 0
            )

        sm = stats(sm_cells)
        vp = stats(vp_cells)
        lactate = self._get_field_array(self.field.lactate)
        self._log(f"\n=== BIOFILM STATUS [MCS={mcs}] ===")
        self._log(f"S. mutans: n={sm['count']}, avgV={sm['avg_volume']:.1f}, totV={sm['total_volume']:.0f}, "
                  f"avgM={sm['avg_mass']:.2e}, totM={sm['total_mass']:.2e}")
        self._log(f"V. parvula: n={vp['count']}, avgV={vp['avg_volume']:.1f}, totV={vp['total_volume']:.0f}, "
                  f"avgM={vp['avg_mass']:.2e}, totM={vp['total_mass']:.2e}")
        self._log(f"Metabolites: lactate max={np.max(lactate):.4f}, avg={np.mean(lactate):.4f}")
        self._log("=" * 40)

    def _setup_plots(self):
        try:
            self.pop_win = self.add_new_plot_window(
                title='Bacterial Populations', x_axis_title='MCS', y_axis_title='Cell Count')
            self.pop_win.add_plot("S_mutans", style='Dots', color='red', size=3)
            self.pop_win.add_plot("V_parvula", style='Dots', color='blue', size=3)
            self.pop_win.add_data_point("S_mutans", 0, 0)
            self.pop_win.add_data_point("V_parvula", 0, 0)

            self.vol_win = self.add_new_plot_window(
                title='Total Cell Volumes', x_axis_title='MCS', y_axis_title='Total Volume (voxels)')
            self.vol_win.add_plot("S_mutans_vol", style='Dots', color='red', size=3)
            self.vol_win.add_plot("V_parvula_vol", style='Dots', color='blue', size=3)
            self.vol_win.add_data_point("S_mutans_vol", 0, 0)
            self.vol_win.add_data_point("V_parvula_vol", 0, 0)

            self.mass_win = self.add_new_plot_window(
                title='Total Bacterial Mass', x_axis_title='MCS', y_axis_title='Total Mass (kg)')
            self.mass_win.add_plot("S_mutans_mass", style='Dots', color='red', size=3)
            self.mass_win.add_plot("V_parvula_mass", style='Dots', color='blue', size=3)
            self.mass_win.add_data_point("S_mutans_mass", 0, 0)
            self.mass_win.add_data_point("V_parvula_mass", 0, 0)

            self.metabolite_win = self.add_new_plot_window(
                title='Average Metabolite Concentrations', x_axis_title='MCS', y_axis_title='Concentration')
            self.metabolite_win.add_plot("Sucrose", style='Dots', color='green', size=3)
            self.metabolite_win.add_plot("Lactate", style='Dots', color='orange', size=3)
            self.metabolite_win.add_data_point("Sucrose", 0, 0)
            self.metabolite_win.add_data_point("Lactate", 0, 0)

            self._log("Plot windows created successfully")

        except Exception as e:
            self._log(f"Error creating plot windows: {e}")
            import traceback
            traceback.print_exc()

    def _update_plots(self, mcs):
        try:
            sm_n = vp_n = 0
            sm_vol = vp_vol = 0.0
            sm_mass = vp_mass = 0.0

            for cell in self.cell_list:
                if cell.type == params.S_MUTANS:
                    sm_n += 1
                    sm_vol += cell.targetVolume
                    sm_mass += cell.dict.get('mass', params.volume_to_mass(cell.targetVolume))
                elif cell.type == params.V_PARVULA:
                    vp_n += 1
                    vp_vol += cell.targetVolume
                    vp_mass += cell.dict.get('mass', params.volume_to_mass(cell.targetVolume))

            self.pop_win.add_data_point("S_mutans", mcs, float(sm_n))
            self.pop_win.add_data_point("V_parvula", mcs, float(vp_n))

            self.vol_win.add_data_point("S_mutans_vol", mcs, float(sm_vol))
            self.vol_win.add_data_point("V_parvula_vol", mcs, float(vp_vol))

            self.mass_win.add_data_point("S_mutans_mass", mcs, float(sm_mass))
            self.mass_win.add_data_point("V_parvula_mass", mcs, float(vp_mass))

            # Simple centerline average over z
            suc = []
            lac = []
            for z in range(params.nz):
                suc.append(float(self.field.sucrose[params.nx//2, params.ny//2, z]))
                lac.append(float(self.field.lactate[params.nx//2, params.ny//2, z]))
            self.metabolite_win.add_data_point("Sucrose", mcs, float(np.mean(suc)))
            self.metabolite_win.add_data_point("Lactate", mcs, float(np.mean(lac)))

            if mcs % (10 * self.bio_frequency) == 0 and mcs > 0:
                self._log(f"Plot MCS {mcs}: SM={sm_n}(V={sm_vol:.0f}), VP={vp_n}(V={vp_vol:.0f}), "
                          f"Suc={np.mean(suc):.3f}, Lac={np.mean(lac):.3f}")

        except Exception as e:
            self._log(f"Error updating plots at MCS {mcs}: {e}")
            import traceback
            traceback.print_exc()

    def finish(self):
        try:
            self._log("\n=== FINAL TIMING STATISTICS ===")
            for key, times in self.stats.items():
                if times:
                    self._log(f"  {key:20s}: Total={sum(times):.2f}s, Avg={np.mean(times):.3f}s")

            out = self.output_dir
            self.pop_win.save_plot_as_png(os.path.join(out, "populations.png"), 1200, 800)
            self.vol_win.save_plot_as_png(os.path.join(out, "volumes.png"), 1200, 800)
            self.mass_win.save_plot_as_png(os.path.join(out, "masses.png"), 1200, 800)
            self.metabolite_win.save_plot_as_png(os.path.join(out, "metabolites.png"), 1200, 800)
            self._log("Plots saved.")
            self._print_biofilm_status_with_volumes(self.mcs)
            self._log("=== FINISH COMPLETE ===")
        except Exception as e:
            self._log(f"Error in finish: {e}")