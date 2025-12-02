import numpy as np
from neuron import h
from typing import Iterable, Callable, Tuple, Optional

HocSec = type(h.Section())


def _inside(a: int, n: int) -> bool:
    return 0 <= a < n


def _idx(i: int, j: int, k: int, Ny: int, Nz: int) -> int:
    return (i * Ny + j) * Nz + k


def _k_rate(diff_um2_ms: float, dx_um: float) -> float:
    """Finite-difference coupling rate (1/ms) for neighbor: D / dx^2."""
    return float(diff_um2_ms) / (float(dx_um) ** 2)

def record_voxels(voxel_indices, voxels):
    """Return (time_vec, dict idx->Vector) to record conc over time."""
    from neuron import h
    tvec = h.Vector(); tvec.record(h._ref_t)
    rec = {}
    for idx in voxel_indices:
        v = h.Vector().record(voxels[idx]._ref_conc)
        rec[idx] = v
    return tvec, rec


def attach_no_voxel_lattice(
        voxel_size=11,    # size of one voxel in um
        lattice_size=110,  # size of the entire 3d lattice (assuming cube)
        d_phys=3.3,       # diffusion constant, D (um^2/ms)
        lam_ms=1000,        # decay constant in ms for simplicity (converted later)
        center_initial=None,  # Value to set initial concentration in center voxel nM
        trapezoidal_ref=True,  # trapezoidal synthesis func from exp reference
        record_line=True  # Record straight line through the cibe (for now its x pos dim)
):
    voxels = {}  # {(x,y,z): hoc no_voxel}
    box = (lattice_size, lattice_size, lattice_size)
    offsets = {
        (voxel_size, 0, 0): 'conc_xp', (-voxel_size, 0, 0): 'conc_xn',
        (0, voxel_size, 0): 'conc_yp', (0, -voxel_size, 0): 'conc_yn',
        (0, 0,  voxel_size): 'conc_zp', (0, 0, -voxel_size): 'conc_zn',
    }  # Create a grid of offsets to wire neighbor pointers for voxels

    # Create a dummy section to host all the voxels
    host_sec = h.Section(name='no_anchor')
    xs = list(range(0, box[0]+1, int(voxel_size)))
    ys = list(range(0, box[1]+1, int(voxel_size)))
    zs = list(range(0, box[2] + 1, int(voxel_size)))
    for x in xs:
        for y in ys:
            for z in zs:
                pp = h.no_voxel(host_sec(0.5))
                voxels[(x,y,z)] = pp

    # link neighbors
    for (x,y,z), v in voxels.items():
        for (dx,dy,dz), pname in offsets.items():
            neighbor_idx = (x+dx, y+dy, z+dz)
            target = voxels.get(neighbor_idx, v)
            h.setpointer(target._ref_conc, pname, v)

    # set diffusion params for NO
    lam = np.log(2)/lam_ms
    d = d_phys/(voxel_size**2)
    for v in voxels.values():
        v.dx_pos = v.dx_neg = v.dy_pos = v.dy_neg = v.dz_pos = v.dz_neg = d
        v.lam = lam

    cx, cy, cz = xs[len(xs) // 2], ys[len(ys) // 2], zs[len(zs) // 2]
    center_key = (cx, cy, cz)

    # Option to set an initial concentration for the center voxel
    if center_initial:
        voxels[center_key].conc0 = center_initial

    if trapezoidal_ref:
        tvec = h.Vector([0, 420, 470, 570, 620, 1000000000000])
        fvec = h.Vector([0,   0, 2.5, 2.5,   0, 0])  # nM/ms
        fvec.play(voxels[center_key]._ref_F, tvec, 1)

    if record_line:
        conc_vecs = {}
        for x in xs:
            if x >= cx:
                idx = (x, cy, cz)
                v = h.Vector().record(voxels[idx]._ref_conc)
                conc_vecs[idx] = v

    return {
        'hSec': host_sec,
        'voxels': voxels,
        'conc_vecs': conc_vecs,
        'f_funciton': fvec,
        'tvec': tvec
            }



def attach_no_voxel_lattice_complex(
    somas: Iterable[HocSec],
    coord_fn: Callable[[HocSec], Tuple[float, float, float]],
    *,
    edge_um: float = 110.0,        # cube edge (µm)
    n_side: int = 11,              # 11 x 11 x 11
    cube_origin: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    # no_voxel specifics (from your MOD)
    pp_class_name: str = "no_voxel",
    conc_range_name: str = "conc",
    # neighbor POINTER field names in your MOD:
    neighbor_fields = ("conc_xp", "conc_xn", "conc_yp", "conc_yn", "conc_zp", "conc_zn"),
    # physics (defaults; change as needed)
    D_um2_ms: float = 3.3,         # diffusion coefficient (µm^2/ms)
    lam_1_ms: float = 0.001,         # decay rate λ (1/ms)
    F_nM_ms: float = 0.0,          # source term F (nM/ms)
    # boundary handling
    reflect_boundaries: bool = False,   # use self-pointer for edges → zero-flux
    # optional: wire a cell-side POINTER to read voxel conc
    cell_ptr_name: Optional[str] = None,   # e.g., 'no_ext' if your cell mech declares POINTER no_ext
    cell_reader_mech: Optional[str] = None # e.g., insert('no_reader') before wiring, if needed
):
    """
    Build and wire an N^3 lattice of `no_voxel` PPs that solve:

        conc' = Σ_axis (k_pos * (nbr_pos - conc) + k_neg * (nbr_neg - conc)) - lam*conc + F

    where k_* are set to D / dx^2 (units 1/ms). Connectivity of the original model is untouched.

    Returns a dict with 'voxels', 'dx', 'n_side', 'origin'.
    """
    Nx = Ny = Nz = int(n_side)
    dx = float(edge_um) / float(n_side)
    ox, oy, oz = cube_origin

    # rates (1/ms) for each axis (isotropic here, but you can split per-axis if needed)
    k = _k_rate(D_um2_ms, dx)

    # anchor section to host PPs
    anchor = h.Section(name="no_anchor"); anchor.L = 1; anchor.diam = 1
    pp_ctor = getattr(h, pp_class_name)

    # allocate voxels
    voxels = [None] * (Nx * Ny * Nz)
    for i in range(Nx):
        for j in range(Ny):
            for k3 in range(Nz):
                pp = pp_ctor(anchor(0.5))
                # set scalar params (can be changed later voxel-by-voxel if needed)
                pp.lam = lam_1_ms
                pp.F   = F_nM_ms
                # initialize symmetric neighbor rates (isotropic)
                pp.dx_pos = k; pp.dx_neg = k
                pp.dy_pos = k; pp.dy_neg = k
                pp.dz_pos = k; pp.dz_neg = k
                voxels[_idx(i, j, k3, Ny, Nz)] = pp

    # neighbor wiring
    cref = f"_ref_{conc_range_name}"
    fn_xp, fn_xn, fn_yp, fn_yn, fn_zp, fn_zn = neighbor_fields

    for i in range(Nx):
        for j in range(Ny):
            for k3 in range(Nz):
                me = voxels[_idx(i, j, k3, Ny, Nz)]

                def link(di, dj, dk, field):
                    ii, jj, kk = i + di, j + dj, k3 + dk
                    if _inside(ii, Nx) and _inside(jj, Ny) and _inside(kk, Nz):
                        src = getattr(voxels[_idx(ii, jj, kk, Ny, Nz)], cref)
                    else:
                        if reflect_boundaries:
                            # self-pointer → (nbr - conc) = 0 → zero-flux boundary
                            src = getattr(me, cref)
                        else:
                            # absorbing: set rate 0 by zeroing the matching coefficient
                            # and point to self to keep pointer valid
                            src = getattr(me, cref)
                            # zero the matching rate on this boundary
                            if   field == "conc_xp": me.dx_pos = 0.0
                            elif field == "conc_xn": me.dx_neg = 0.0
                            elif field == "conc_yp": me.dy_pos = 0.0
                            elif field == "conc_yn": me.dy_neg = 0.0
                            elif field == "conc_zp": me.dz_pos = 0.0
                            elif field == "conc_zn": me.dz_neg = 0.0
                    h.setpointer(src, field, me)

                link(+1,  0,  0, fn_xp)
                link(-1,  0,  0, fn_xn)
                link( 0, +1,  0, fn_yp)
                link( 0, -1,  0, fn_yn)
                link( 0,  0, +1, fn_zp)
                link( 0,  0, -1, fn_zn)

    # helper to map coords → voxel index
    def to_index(x, y, z):
        ix = int((x - ox) // dx); iy = int((y - oy) // dx); iz = int((z - oz) // dx)
        ix = max(0, min(Nx - 1, ix)); iy = max(0, min(Ny - 1, iy)); iz = max(0, min(Nz - 1, iz))
        return ix, iy, iz

    # wire each soma to its voxel concentration, if requested
    if cell_ptr_name is not None:
        for sec in somas:
            if cell_reader_mech:
                try: sec.insert(cell_reader_mech)
                except: pass
            x, y, z = coord_fn(sec)
            ix, iy, iz = to_index(x, y, z)
            vpp = voxels[_idx(ix, iy, iz, Ny, Nz)]
            h.setpointer(getattr(vpp, cref), cell_ptr_name, sec(0.5))

    return {
        "voxels": voxels,
        "dx": dx,
        "n_side": n_side,
        "origin": (ox, oy, oz),
    }
