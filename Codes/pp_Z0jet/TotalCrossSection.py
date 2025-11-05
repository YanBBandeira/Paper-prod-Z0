import numpy as np
import vegas
from scipy.interpolate import RegularGridInterpolator
import os

# ======================================================
# 1. Physical constants and decay model
# ======================================================
class PhysicsParameters:
    def __init__(self):
        self.pi = np.pi
        self.alfem = 1/137.0
        self.sin2 = 0.23
        self.aw = np.arcsin(np.sqrt(self.sin2))
        self.Mz = 91.2  # GeV
        self.rs = 13000.0  # sqrt(s)
        self.ml = 0.1056  # muon mass [GeV]

    def dilepton_decay(self, M):
        """Partial width × branching ratio for Z → μ⁺μ⁻."""
        Mz2 = self.Mz**2
        M2 = M**2
        width = ((self.alfem*M)/(6.0*(np.sin(2*self.aw))**2)) * (
            (160.0/3.0)*(np.sin(self.aw)**4) - 40.0*(np.sin(self.aw)**2) + 21.0
        )
        branch = 3.3 / 100.0
        inv_mass_dist = (1.0/np.pi) * ((M*width)/((M2 - Mz2)**2 + (M*width)**2))
        return inv_mass_dist * branch


# ======================================================
# 2. Histogram manager
# ======================================================
class HistogramManager:
    def __init__(self, y_bins, pt_bins, m_bins):
        self.y_edges = np.linspace(*y_bins)
        self.pt_edges = np.linspace(*pt_bins)
        self.m_edges = np.linspace(*m_bins)

        self.sig_y = np.zeros(len(self.y_edges)-1)
        self.sig_pt = np.zeros(len(self.pt_edges)-1)
        self.sig_m = np.zeros(len(self.m_edges)-1)

        self.y_slices = [(2.0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.0), (4.0, 4.5)]
        self.sig_ypt = [np.zeros(len(self.pt_edges)-1) for _ in self.y_slices]

    def fill(self, y, pt, M, weight):
        # y distribution
        iy = np.digitize(y, self.y_edges) - 1
        if 0 <= iy < len(self.sig_y):
            self.sig_y[iy] += weight

        # pT distribution
        ipt = np.digitize(pt, self.pt_edges) - 1
        if 0 <= ipt < len(self.sig_pt):
            self.sig_pt[ipt] += weight

        # M distribution
        im = np.digitize(M, self.m_edges) - 1
        if 0 <= im < len(self.sig_m):
            self.sig_m[im] += weight

        # y-pT double differential
        for idx, (ymin, ymax) in enumerate(self.y_slices):
            if ymin < y < ymax and 0 <= ipt < len(self.pt_edges)-1:
                deltaY = ymax - ymin
                self.sig_ypt[idx][ipt] += weight / deltaY

    def write_results(self, outdir="Output"):
        os.makedirs(outdir, exist_ok=True)
        np.savetxt(f"{outdir}/dsig_dy.dat",
                   np.column_stack([0.5*(self.y_edges[1:]+self.y_edges[:-1]), self.sig_y]))
        np.savetxt(f"{outdir}/dsig_dpt.dat",
                   np.column_stack([0.5*(self.pt_edges[1:]+self.pt_edges[:-1]), self.sig_pt]))
        np.savetxt(f"{outdir}/dsig_dm.dat",
                   np.column_stack([0.5*(self.m_edges[1:]+self.m_edges[:-1]), self.sig_m]))
        for i, (ymin, ymax) in enumerate(self.y_slices):
            np.savetxt(f"{outdir}/dsig_dydpt_y{ymin:.1f}-{ymax:.1f}.dat",
                       np.column_stack([0.5*(self.pt_edges[1:]+self.pt_edges[:-1]),
                                        self.sig_ypt[i]]))
        print(f"Histograms written to {outdir}/")


# ======================================================
# 3. Grid interpolator
# ======================================================
class GridInterpolator:
    def __init__(self, filename, n_points=15):
        self.n_points = n_points
        self.filename = filename
        self.y_grid = np.zeros(n_points)
        self.pt_grid = np.zeros(n_points)
        self.m_grid = np.zeros(n_points)
        self.parton_grid = np.zeros((n_points, n_points, n_points))
        self._read_grid()

    def _read_grid(self):
        with open(self.filename, "r") as f:
            for i in range(self.n_points):
                for j in range(self.n_points):
                    for k in range(self.n_points):
                        line = f.readline()
                        if not line:
                            break
                        yv, ptv, mv, val = map(float, line.strip().split())
                        self.y_grid[i] = yv
                        self.pt_grid[j] = ptv
                        self.m_grid[k] = mv
                        self.parton_grid[i, j, k] = val

    def interpolate(self, y, pt, m):
        interpolator = RegularGridInterpolator(
            (self.y_grid, self.pt_grid, self.m_grid),
            self.parton_grid,
            bounds_error=False,
            fill_value=None
        )
        return float(interpolator((y, pt, m)))


# ======================================================
# 4. Cross section integrand with full kinematics
# ======================================================

# ======================================================
# 4. Cross section integrand with full Fortran-style kinematics
# ======================================================
class CrossSectionIntegrand:
    def __init__(self, params, hist_manager, hadronic):
        self.params = params
        self.hist = hist_manager
        self.hadronic = hadronic

    def __call__(self, x):
        yp, ym, ktp, ktm, phip, phim = x
        mp = mm = self.params.ml

        # === CORTES ===
        if not (2.0 <= yp <= 4.5 or 2.0 <= ym <= 4.5):
            return 0.0
        if not (20.0 <= ktp or 20.0 <= ktm):
            return 0.0
        # -----------------------------------------------
        # Componentes transversas dos múons
        # -----------------------------------------------
        ktpx = ktp * np.cos(phip)
        ktpy = ktp * np.sin(phip)
        ktmx = ktm * np.cos(phim)
        ktmy = ktm * np.sin(phim)

        # -----------------------------------------------
        # Massas transversas dos múons (m⊥ = √(kT² + mμ²))
        # -----------------------------------------------
        mperp_p2 = ktp**2 + mp**2
        mperp_m2 = ktm**2 + mm**2
        mperp_p = np.sqrt(mperp_p2)
        mperp_m = np.sqrt(mperp_m2)

        # -----------------------------------------------
        # Diferença de rapidez e momento transversal do bóson
        #    (seguindo o sinal em pty do código Fortran)
        # -----------------------------------------------
        ptx = ktpx + ktmx
        pty = ktpy - ktmy   # <- sinal trocado conforme Fortran
        pt2 = ptx**2 + pty**2
        pt = np.sqrt(pt2)

        # -----------------------------------------------
        # Massa invariante do sistema (Fortran-style)
        #     m² = m⊥₊² + m⊥₋² + 2·m⊥₊·m⊥₋·cosh(Δy) - pT²
        # -----------------------------------------------
        deltaY = yp - ym
        m2 = mperp_p**2 + mperp_m**2 + 2.0 * mperp_p * mperp_m * np.cosh(deltaY) - pt2
        if m2 <= 0:
            return 0.0
        M = np.sqrt(m2)

        # -----------------------------------------------
        # Rapidez do bóson (Fortran logic)
        #     y = log( (xf * rs) / sqrt(pt² + M²) )
        #     onde xf = xp + xm
        # -----------------------------------------------
        rs = self.params.rs
        xp = (ktp / rs) * np.exp(yp)
        xm = (ktm / rs) * np.exp(ym)
        xf = xp + xm
        y = np.log(xf * (rs / np.sqrt(pt2 + M**2)))

        # -----------------------------------------------
        # Frações de momento partônico (Fortran form)
        #     x1 = √(M² + pT²)/rs * e^{+y}
        #     x2 = √(M² + pT²)/rs * e^{−y}
        # -----------------------------------------------
        x1 = np.sqrt(M**2 + pt**2) / rs * np.exp(+y)
        x2 = np.sqrt(M**2 + pt**2) / rs * np.exp(-y)

        # -----------------------------------------------
        # Cortes cinemáticos (como no Fortran)
        # -----------------------------------------------
        if not (2.0 <= y <= 4.5):
            return 0.0
        if not (60.0 <= M <= 120.0):
            return 0.0
        if x1 >= 1.0 or x2 >= 1.0:
            return 0.0

        # -----------------------------------------------
        # Jacobiano e pre-integral (Fortran-style)
        # -----------------------------------------------
        varJacobian = (2.0 / rs) * np.sqrt(M**2 + pt2) * np.cosh(y)
        preIntegral = (x1 / (x1 + x2)) * varJacobian
        
        
        # -----------------------------------------------
        # Fatores hadrônicos e decaimento do bóson
        # -----------------------------------------------
        hadronic_val = self.hadronic.interpolate(y, pt, M)
        decay = self.params.dilepton_decay(M)

        # Integrando o pre-integral
        sigma = preIntegral * decay * hadronic_val

        # -----------------------------------------------
        # Preenchimento dos histogramas com o peso do VEGAS
        # -----------------------------------------------
        self.hist.fill(y, pt, M, sigma)
        return sigma


    

# ======================================================
# 5. VEGAS integrator
# ======================================================
class VegasIntegrator:
    def __init__(self, integrand, bounds):
        self.integrand = integrand
        self.bounds = bounds
        self.integrator = vegas.Integrator(bounds)

    def integrate(self, nitn=10, neval=10000):
        result = self.integrator(self.integrand, nitn=nitn, neval=neval)
        print(f"VEGAS ⟨σ⟩ = {result.mean:.6e} ± {result.sdev:.6e}  (Q={result.Q:.2f})")
        return result.mean, result.sdev


# ======================================================
# 6. Main driver
# ======================================================
def main():
    params = PhysicsParameters()

    hist = HistogramManager(
        y_bins=(2.0, 4.5, 61),
        pt_bins=(0.0, 150.0, 61),
        m_bins=(60.0, 120.0, 61)
    )

    hadronic = GridInterpolator("hadronic_grid.dat", n_points=15)

    bounds = [
        (2.0, 4.5),       # y+
        (2.0, 4.5),       # y-
        (20.0, 120.0),    # kT+
        (20.0, 120.0),    # kT-
        (-np.pi, np.pi),  # phi+
        (-np.pi, np.pi)   # phi-
    ]

    integrand = CrossSectionIntegrand(params, hist, hadronic)
    vegas_integrator = VegasIntegrator(integrand, bounds)

    # Warm-up
    vegas_integrator.integrate(nitn=5, neval=5000)
    # Production run
    vegas_integrator.integrate(nitn=10, neval=20000)

    hist.write_results()


if __name__ == "__main__":
    main()
