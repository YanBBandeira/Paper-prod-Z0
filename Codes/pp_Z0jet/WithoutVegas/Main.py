import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os

#1. Physical constants and decay model
class PhysicsParameters:
    def __init__(self):
        self.pi = np.pi 
        self.alfem = 1/137.0
        self.sin2 = 0.23
        self.aw = np.arcsin(np.sqrt(self.sin2))
        self.Mz = 91.2 #GeV
        self.rs = 13000.0 # sqrt(s)
        self.ml = 0.1056 # muon mass [GeV]
        
    def dilepton_decay(self, M):
        Mz2 = self.Mz**2.0
        M2 = M**2
        width = ((self.alfem*M)/(6.0*(np.sin(2*self.aw)**2.0))) * (
            (160.0/3.0)*(np.sin(self.aw)**4.0) - 40.0*(np.sin(self.aw)**2.0) + 21.0
        )
        branch = 3.3 / 100.0
        inv_mass_dist = (1.0/np.pi) * (
            (M * width)/((M2 - Mz2)**2.0 + (M * width)**2.0)
        )
        return inv_mass_dist * branch
    
#2. Histogram manager
class HistogramManager:
    def __init__(self, y_bins, pt_range, m_bins):
        self.y_edges = np.linspace(*y_bins)
        self.pt_edges = np.logspace(np.log10(pt_range[0]), np.log10(pt_range[1]), pt_range[2])
        self.m_edges = np.linspace(*m_bins)
        
        self.sig_y = np.zeros(len(self.y_edges)-1)
        self.sig_pt = np.zeros(len(self.pt_edges)-1)
        self.sig_m = np.zeros(len(self.m_edges)-1)
        
        self.y_slices = [(2.0,2.5),(2.5,3.0),(3.0,3.5),(3.5,4.0),(4.0,4.5)]
        self.sig_ypt = [np.zeros(len(self.pt_edges)-1) for _ in self.y_slices] # List with the sliced ditros
        
    def fill(self, y,pt,M,weight):
        # y distribution
        iy = np.digitize(y, self.y_edges) - 1
        if 0 <= iy < len(self.sig_y):
            self.sig_y[iy] += weight
            
        # pT distribution
        ipt = np.digitize(pt, self.pt_edges) - 1
        if  0<= ipt < len(self.sig_pt):
            self.sig_pt[ipt] += weight
            
        # M distribution
        im = np.digitize(M, self.m_edges) - 1
        if 0 <= im < len(self.sig_m):
            self.sig_m[im] += weight
            
        # y-pT double differential
        for idx, (ymin, ymax) in enumerate(self.y_slices):
            if ymin < y < ymax and 0 <= ipt < len(self.pt_edges) -1:
                deltaY = ymax - ymin
                self.sig_ypt[idx][ipt] += weight / deltaY
    
    def write_results(self, outdir=r"C:\Users\Callidus\Documents\Clones\Paper-prod-Z0\Codes\pp_Z0jet\Output"):
        os.makedirs(outdir, exist_ok=True)
        
        dy = np.diff(self.y_edges)
        dpt = np.diff(self.pt_edges)
        dm = np.diff(self.m_edges)
        
        np.savetxt(f"{outdir}/dsig_dy.dat",
                    np.column_stack([
                        0.5*(self.y_edges[1:] + self.y_edges[:-1]), self.sig_y / dy
                                    ])
                   )  
        
        np.savetxt(f"{outdir}/dsig_dpt.dat",
                    np.column_stack([
                        0.5*(self.pt_edges[1:] + self.pt_edges[:-1]), self.sig_pt / dpt
                                    ])
                   )  
        
        np.savetxt(f"{outdir}/dsig_dm.dat",
                    np.column_stack([
                        0.5*(self.m_edges[1:] + self.m_edges[:-1]), self.sig_m / dm
                                    ])
                   )            
        
        for i, (ymin, ymax) in enumerate(self.y_slices):
            np.savetxt(f"{outdir}/dsig_dydpt_y{ymin: .0f}-{ymax: .0f}.dat",
                       np.column_stack([
                           0.5*(self.pt_edges[1:] + self.pt_edges[:-1]),
                           self.sig_ypt[i] / dpt
                       ])
                       )
        
        print(f"Histograms written to {outdir}/")
        
 # 3. Grid interpolator
class GridInterpolator:
    def __init__(self, filename, n_points):
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
    
# 4. Cross section integrand with full kinematics
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
             return 0.0, 0.0, 0.0, 0.0
        if not (20.0 <= ktp or 20.0 <= ktm):
            return 0.0, 0.0, 0.0, 0.0

        ktpx = ktp * np.cos(phip)
        ktpy = ktp * np.sin(phip)
        ktmx = ktm * np.cos(phim)
        ktmy = ktm * np.sin(phim)

        mperp_p2 = ktp**2 + mp**2
        mperp_m2 = ktm**2 + mm**2
        mperp_p = np.sqrt(mperp_p2)
        mperp_m = np.sqrt(mperp_m2)

        ptx = ktpx + ktmx
        pty = ktpy - ktmy
        pt2 = ptx**2 + pty**2
        pt = np.sqrt(pt2)

        deltaY = yp - ym
        m2 = mperp_p**2 + mperp_m**2 + 2.0 * mperp_p * mperp_m * np.cosh(deltaY) - pt2
        if m2 <= 0:
             return 0.0, 0.0, 0.0, 0.0
        M = np.sqrt(m2)

        rs = self.params.rs
        xp = (ktp / rs) * np.exp(yp)
        xm = (ktm / rs) * np.exp(ym)
        xf = xp + xm
        y = np.log(xf * (rs / np.sqrt(pt2 + M**2)))

        x1 = np.sqrt(M**2 + pt**2) / rs * np.exp(+y)
        x2 = np.sqrt(M**2 + pt**2) / rs * np.exp(-y)

        if not (60.0 <= M <= 120.0):
            return 0.0, 0.0, 0.0, 0.0
        if x1 >= 1.0 or x2 >= 1.0:
             return 0.0, 0.0, 0.0, 0.0

        varJacobian = (2.0 / rs) * np.sqrt(M**2 + pt2) * np.cosh(y)
        preIntegral = (x1 / (x1 + x2)) * varJacobian

        hadronic_val = self.hadronic.interpolate(y, pt, M)
        decay = self.params.dilepton_decay(M)

        sigma = decay * hadronic_val * preIntegral * ktm * ktp

        #self.hist.fill(y, pt, M, sigma)
        #print(f"Valor do sigma {sigma}, do pt {pt} e da rapidez {y}")
        return sigma, y, pt, M
        
# 5. Grids integration
class GridIntegrator:
    def __init__(self, integrand, hist_manager, nPoints=20,
                 phiP_limits=(0.0, np.pi),
                 phiM_limits=(0.0, np.pi),
                 kP_limits=(20.0, 200.0),
                 kM_limits=(20.0, 200.0),
                 yP_limits=(2.0, 4.5),
                 yM_limits=(2.0, 4.5),
                 ):
        self.integrand = integrand
        self.nPoints = nPoints
        self.hist = hist_manager

        # Recebe os limites como tuplas (min, max)
        self.phiPmin, self.phiPmax = phiP_limits
        self.phiMmin, self.phiMmax = phiM_limits
        self.kPmin, self.kPmax = kP_limits
        self.kMmin, self.kMmax = kM_limits
        self.yPmin, self.yPmax = yP_limits
        self.yMmin, self.yMmax = yM_limits

        # Calcular os passos
        self.dphiP = (self.phiPmax - self.phiPmin) / self.nPoints
        self.dphiM = (self.phiMmax - self.phiMmin) / self.nPoints
        self.dkP = (self.kPmax - self.kPmin) / self.nPoints
        self.dkM = (self.kMmax - self.kMmin) / self.nPoints
        self.dyP = (self.yPmax - self.yPmin) / self.nPoints
        self.dyM = (self.yMmax - self.yMmin) / self.nPoints

    def integrate(self):
        total_integral = 0.0

        for iphiP in range(self.nPoints):
            phiP = self.phiPmin + (iphiP + 0.5) * self.dphiP
            for iphiM in range(self.nPoints):
                phiM = self.phiMmin + (iphiM + 0.5) * self.dphiM
                for ikP in range(self.nPoints):
                    kP = self.kPmin + (ikP + 0.5) * self.dkP
                    for ikM in range(self.nPoints):
                        kM = self.kMmin + (ikM + 0.5) * self.dkM
                        for iyP in range(self.nPoints):
                            yP = self.yPmin + (iyP + 0.5) * self.dyP
                            for iyM in range(self.nPoints):
                                yM = self.yMmin + (iyM + 0.5) * self.dyM

                                x = (yP, yM, kP, kM, phiP, phiM)
                                sigma, y, pt, M = self.integrand(x)

                                dV = (self.dphiP * self.dphiM * self.dkP *
                                      self.dkM * self.dyP * self.dyM)

                                total_integral += sigma * dV
                                self.integrand.hist.fill(y, pt, M, sigma * dV)

        print(f"Integral estimada pelo grid: {total_integral:.6e}")
        return total_integral 
   
   
# 6. Main
def main():
    params = PhysicsParameters()

    hist = HistogramManager(
        y_bins=(2.0, 4.5, 61),
        pt_range=(1.0, 150.0, 61),
        m_bins=(60.0, 120.0, 61)
    )

    hadronic = GridInterpolator(r"C:\Users\Callidus\Documents\Clones\Paper-prod-Z0\Codes\pp_Z0jet\Grids\DatFiles\tst_grid.dat", n_points=25)

    integrand = CrossSectionIntegrand(params, hist, hadronic)

    grid_integrator = GridIntegrator(
        integrand,
        hist,
        nPoints=10,
        phiP_limits=(0.0, np.pi),
        phiM_limits=(0.0, np.pi),
        kP_limits=(20.0, 200.0),
        kM_limits=(20.0, 200.0),
        yP_limits=(2.0, 4.5),
        yM_limits=(2.0, 4.5)     
    )
    grid_integrator.integrate()

    hist.write_results()
    
if __name__ == "__main__":
    main()