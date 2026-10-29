"""
OpenMM Molecular Dynamics engine: energy minimization, NVT thermal equilibration,
and Steered Molecular Dynamics (SMD) mechanical pulling simulation.
"""

from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass
import io
import numpy as np
import openmm as mm
from openmm import app, unit


# Conversion constant: 1 kJ/(mol * nm) = 1.66054 pN
KJ_MOL_NM_TO_PN = 1.660539


@dataclass
class SMDFrame:
    step: int
    time_ps: float
    end_to_end_nm: float
    target_pos_nm: float
    pulled_atom_pos_nm: float
    force_pn: float
    work_kj_mol: float
    potential_energy_kj_mol: float


@dataclass
class SMDResult:
    sequence: str
    num_residues: int
    total_steps: int
    total_time_ps: float
    pull_velocity_nm_ps: float
    spring_k_kj_mol_nm2: float
    frames: List[SMDFrame]
    initial_distance_nm: float
    final_distance_nm: float
    peak_force_pn: float
    total_work_kj_mol: float
    trajectory_pdb: str


def setup_simulation_system(
    pdb_text: str,
    forcefield_names: Tuple[str, ...] = ("amber99sbildn.xml", "amber99_obc.xml"),
    temperature_k: float = 300.0,
    friction_per_ps: float = 1.0,
    time_step_fs: float = 2.0,
    platform_name: str = "CPU",
) -> Tuple[app.Simulation, app.PDBFile, List[int]]:
    """
    Set up OpenMM simulation with implicit solvent and identify CA atoms.
    """
    pdb = app.PDBFile(io.StringIO(pdb_text))
    forcefield = app.ForceField(*forcefield_names)
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )

    ca_indices = [
        atom.index for atom in pdb.topology.atoms() if atom.name == "CA"
    ]
    if len(ca_indices) < 2:
        raise ValueError("System must have at least 2 alpha carbons (CA).")

    integrator = mm.LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        friction_per_ps / unit.picosecond,
        time_step_fs * unit.femtoseconds,
    )

    # Use CPU or Reference platform
    try:
        platform = mm.Platform.getPlatformByName(platform_name)
    except Exception:
        platform = mm.Platform.getPlatformByName("Reference")

    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    return simulation, pdb, ca_indices


def minimize_energy(simulation: app.Simulation, max_iterations: int = 150) -> float:
    """Run conjugate gradient energy minimization; return potential energy in kJ/mol."""
    simulation.minimizeEnergy(maxIterations=max_iterations)
    state = simulation.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)


def run_steered_md(
    pdb_text: str,
    sequence: str,
    steps: int = 1000,
    report_interval: int = 20,
    pull_velocity_nm_ps: float = 0.05,
    spring_k_kj_mol_nm2: float = 500.0,
    temperature_k: float = 300.0,
    minimize_first: bool = True,
    platform_name: str = "CPU",
) -> SMDResult:
    """
    Perform Steered Molecular Dynamics (SMD) mechanical pulling simulation:
    - Anchors N-terminal CA atom
    - Pulls C-terminal CA atom with harmonic spring along X-axis
    - Records tension force, extension, and mechanical work over time
    """
    if steps <= 0:
        raise ValueError("steps must be a positive integer.")
    if report_interval <= 0:
        raise ValueError("report_interval must be a positive integer.")
    if steps % report_interval != 0:
        raise ValueError("steps must be divisible by report_interval so all steps are simulated.")
    if pull_velocity_nm_ps <= 0 or spring_k_kj_mol_nm2 <= 0:
        raise ValueError("pull velocity and spring constant must be positive.")

    simulation, pdb, ca_indices = setup_simulation_system(
        pdb_text=pdb_text,
        temperature_k=temperature_k,
        platform_name=platform_name,
    )

    start_ca = ca_indices[0]
    end_ca = ca_indices[-1]

    # Energy minimization
    if minimize_first:
        minimize_energy(simulation, max_iterations=150)

    initial_pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    start_pos = initial_pos[start_ca]
    end_pos = initial_pos[end_ca]
    d_initial = float(np.linalg.norm(end_pos - start_pos))

    # Add harmonic restraint to anchor N-terminal CA
    anchor_force = mm.CustomExternalForce("1000.0 * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    anchor_force.addPerParticleParameter("x0")
    anchor_force.addPerParticleParameter("y0")
    anchor_force.addPerParticleParameter("z0")
    anchor_force.addParticle(start_ca, [start_pos[0], start_pos[1], start_pos[2]])

    # Add moving harmonic spring to pull C-terminal CA along X axis
    pull_force = mm.CustomExternalForce("0.5 * k_pull * ((x - (x0 + v_pull*t))^2 + (y-y0)^2 + (z-z0)^2)")
    pull_force.addGlobalParameter("k_pull", spring_k_kj_mol_nm2 * unit.kilojoules_per_mole / unit.nanometer**2)
    pull_force.addGlobalParameter("v_pull", pull_velocity_nm_ps * unit.nanometer / unit.picosecond)
    pull_force.addGlobalParameter("t", 0.0 * unit.picosecond)
    pull_force.addPerParticleParameter("x0")
    pull_force.addPerParticleParameter("y0")
    pull_force.addPerParticleParameter("z0")
    pull_force.addParticle(end_ca, [end_pos[0], end_pos[1], end_pos[2]])

    system = simulation.system
    system.addForce(anchor_force)
    system.addForce(pull_force)

    # Re-initialize context to pick up new forces
    current_positions = simulation.context.getState(getPositions=True).getPositions()
    simulation.context.reinitialize(preserveState=True)
    simulation.context.setPositions(current_positions)

    frames: List[SMDFrame] = []
    traj_io = io.StringIO()
    dt_ps = 0.002  # 2 fs timestep
    num_cycles = steps // report_interval
    accumulated_work = 0.0
    prev_x = end_pos[0]
    peak_force = 0.0

    for cycle in range(num_cycles):
        # Update the moving restraint every integration step. Updating only once
        # per report interval would create a staircase pulling protocol.
        for step_offset in range(report_interval):
            current_step = cycle * report_interval + step_offset + 1
            current_time_ps = current_step * dt_ps
            simulation.context.setParameter("t", current_time_ps)
            simulation.step(1)

        current_step = (cycle + 1) * report_interval
        current_time_ps = current_step * dt_ps
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        pot_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

        # Append to trajectory PDB string
        app.PDBFile.writeModel(pdb.topology, state.getPositions(), traj_io, modelIndex=cycle + 1)

        d_current = float(np.linalg.norm(pos[end_ca] - pos[start_ca]))
        target_x = end_pos[0] + (pull_velocity_nm_ps * current_time_ps)
        pulled_x = pos[end_ca][0]

        # Force = k * (target - actual)
        force_kj_mol_nm = spring_k_kj_mol_nm2 * (target_x - pulled_x)
        force_pn = force_kj_mol_nm * KJ_MOL_NM_TO_PN
        peak_force = max(peak_force, force_pn)

        # Work performed by the moving harmonic restraint is F * d(lambda),
        # where lambda is the spring target position (not the atom displacement).
        target_dx = target_x - (end_pos[0] + pull_velocity_nm_ps * ((current_step - report_interval) * dt_ps))
        work_increment = force_kj_mol_nm * target_dx
        accumulated_work += work_increment

        frames.append(
            SMDFrame(
                step=current_step,
                time_ps=round(current_time_ps, 3),
                end_to_end_nm=round(d_current, 4),
                target_pos_nm=round(target_x, 4),
                pulled_atom_pos_nm=round(pulled_x, 4),
                force_pn=round(force_pn, 2),
                work_kj_mol=round(accumulated_work, 3),
                potential_energy_kj_mol=round(pot_energy, 2),
            )
        )

    d_final = frames[-1].end_to_end_nm if frames else d_initial
    total_time = (steps * dt_ps)

    return SMDResult(
        sequence=sequence,
        num_residues=len(sequence),
        total_steps=steps,
        total_time_ps=round(total_time, 3),
        pull_velocity_nm_ps=pull_velocity_nm_ps,
        spring_k_kj_mol_nm2=spring_k_kj_mol_nm2,
        frames=frames,
        initial_distance_nm=round(d_initial, 3),
        final_distance_nm=round(d_final, 3),
        peak_force_pn=round(peak_force, 2),
        total_work_kj_mol=round(accumulated_work, 3),
        trajectory_pdb=traj_io.getvalue(),
    )
