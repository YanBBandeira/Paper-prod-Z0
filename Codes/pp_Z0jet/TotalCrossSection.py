import vegas 
import numpy as np
from scipy.interpolate import RegularGridInterpolator



def main():
    global hist, bin_params
    # Binning parameters
    ny = 100
    y_min = 2.0
    y_max = 4.5
    dy = (y_max - y_min) / ny

    npt = 100
    pt_min = 0.0
    pt_max = 150.0
    dpt = (pt_max - pt_min) / npt

    nm = 100
    m_min = 60.0
    m_max = 120.0
    dm = (m_max - m_min) / nm

    bin_params = {
        'ny': ny, 'y_min': y_min, 'y_max': y_max, 'dy': dy,
        'npt': npt, 'pt_min': pt_min, 'pt_max': pt_max, 'dpt': dpt,
        'nm': nm, 'm_min': m_min, 'm_max': m_max, 'dm': dm
    }

    # Histograms
    hist = {
        'iDoHist': 0,
        'sig_y': np.zeros(ny),
        'sig_pt': np.zeros(npt),
        'sig_m': np.zeros(nm),
        'sig_ypt1': np.zeros(npt),
        'sig_ypt2': np.zeros(npt),
        'sig_ypt3': np.zeros(npt),
        'sig_ypt4': np.zeros(npt),
        'sig_ypt5': np.zeros(npt)
    }
    
    integ = vegas.Integrator([[0, 1], [0, 1]])
    hist['iDoHist'] == 0
    integ(IntegrandSigma, nitn=10, neval=10000)
    
    hist['iDoHist'] == 1
    result = integ(IntegrandSigma, nitn=10, neval=10000)
    print('Resultado:', result.mean, '+-', result.sdev)
    
    # Fill histograms with iDoHist=1
    # Output results
    with open('dsig_dy.dat', 'w') as f:
        sum_y = 0.0
        for iy in range(ny):
            y = y_min + iy * dy - dy / 2.0
            f.write(f"{y:8.4f} {hist['sig_y'][iy]:12.4e}\n")
            sum_y += hist['sig_y'][iy] * dy
        print(f"sum_y: {sum_y:.6e}")

    with open('dsig_dpt.dat', 'w') as f:
        sum_pt = 0.0
        for ipt in range(npt):
            pt = pt_min + ipt * dpt - dpt / 2.0
            f.write(f"{pt:8.4f} {hist['sig_pt'][ipt]:12.4e}\n")
            sum_pt += hist['sig_pt'][ipt] * dpt
        print(f"sum_pt: {sum_pt:.6e}")

    with open('dsig_dm.dat', 'w') as f:
        for im in range(nm):
            m = m_min + im * dm - dm / 2.0
            f.write(f"{m:8.4f} {hist['sig_m'][im]:12.4e}\n")
    
    
if __name__ == "__main__":
    main()


# ===========================
# Parameters (module)
# ===========================
pi = 4.0 * np.arctan(1.0)
pi2 = pi ** 2
alfem = 1.0 / 137.0
sin2 = 0.23
aw = np.arcsin(np.sqrt(sin2))
Mz = 91.2  # Z0 mass
rs = 13000.0  # sqrt(s)

# ===========================
# Dilepton Decay Function
# ===========================
def DileptonDecay(MVar):
    M = MVar
    M2 = M * M
    Mz2 = Mz * Mz
    DecayWidth = ((alfem * M) / (6.0 * (np.sin(2.0 * aw) ** 2))) * (
        (160.0 / 3.0) * (np.sin(aw) ** 4) - 40.0 * (np.sin(aw) ** 2) + 21.0
    )
    Branch = 3.3662 / 100.0
    InvariantMassDist = (1.0 / pi) * (
        (M * DecayWidth) / ((M2 - Mz2) ** 2 + (M * DecayWidth) ** 2)
    )
    Result = InvariantMassDist * Branch
    return Result

def read_grid(nPoints, filename):
    # Inicializa arrays
    yGrid = np.zeros(nPoints)
    ptGrid = np.zeros(nPoints)
    PartonLevelGrid = np.zeros((nPoints, nPoints))

    with open(filename, 'r') as f:
        lines = f.readlines()

    idx = 0
    for i in range(nPoints):
        for j in range(nPoints):
            yVar, ptVar, PartonLevelVar = map(float, lines[idx].split())
            yGrid[i] = yVar
            ptGrid[j] = ptVar
            PartonLevelGrid[i, j] = PartonLevelVar
            idx += 1

    return yGrid, ptGrid, PartonLevelGrid

def InterpolateGrid(y, pt, filename='kslinear_grid.dat', nPoints=5000):
    yGrid, ptGrid, PartonLevelGrid = read_grid(nPoints, filename)
    # Cria interpolador
    interpolator = RegularGridInterpolator((yGrid, ptGrid), PartonLevelGrid)
    value = interpolator([[y, pt]])[0]
    return value

# ===========================
# IntegrandSigma Subroutine
# ===========================
def IntegrandSigma(x):

    pt2Var_max = 150**2
    pt2Var_min = 1e-2**2
    m2Var_min = 60**2
    m2Var_max = 120**2
    yVar_max = 4.5
    yVar_min = 2.0
    yVar     = yVar_min + (yVar_max - yVar_min)*x[0]
    m2Var    = m2Var_min + (m2Var_max - m2Var_min)*x[1]
    pt2Var   = pt2Var_min + (pt2Var_max - pt2Var_min)*x[2]
    physicalWgt = vegas.RAvg*jac/vegas.itermax
    
    global x2, pt, x1, M
    M2 = m2Var
    M  = np.sqrt(m2Var)
    pt = np.sqrt(pt2Var)
    y = yVar

    x1 = (np.sqrt(M2 + pt ** 2) / rs) * np.exp(y)
    x2 = (np.sqrt(M2 + pt ** 2) / rs) * np.exp(-y)

    varJacobian = (2.0 / rs) * np.sqrt(M2 + pt ** 2) * np.cosh(yVar)
    preIntegral = (x1 / (x1 + x2)) * varJacobian

    # HadronicCrossSection = DGAUSS(IntegrandHadronicCrossSection, 0.0, 1.0, 1e-4)
    # Result = preIntegral * DileptonDecay(M) * HadronicCrossSection
    # For demonstration, set HadronicCrossSection = 1.0
    HadronicCrossSection = InterpolateGrid(y, pt)
    Result = preIntegral * DileptonDecay(M) * HadronicCrossSection

    units = 0.389e9  # GeV^-2 to pb
    SigTot = pi * Result * units

    # Boundaries
    if x1 > 1.0 or x2 > 1.0 or M < 60.0 or M > 120.0:
        SigTot = 0.0

    # Indices
    ny, y_min, y_max, dy = bin_params['ny'], bin_params['y_min'], bin_params['y_max'], bin_params['dy']
    npt, pt_min, pt_max, dpt = bin_params['npt'], bin_params['pt_min'], bin_params['pt_max'], bin_params['dpt']
    nm, m_min, m_max, dm = bin_params['nm'], bin_params['m_min'], bin_params['m_max'], bin_params['dm']

    iy = int((y - y_min) / dy)
    ipt = int((pt - pt_min) / dpt)
    im = int((M - m_min) / dm)

    # Collecting spectra
    if hist['iDoHist'] == 1:
        if y_min < y < y_max and 0 <= iy < ny:
            hist['sig_y'][iy] += SigTot * physicalWgt / dy
        if pt_min < pt < pt_max and 0 <= ipt < npt:
            hist['sig_pt'][ipt] += SigTot * physicalWgt / dpt
        if pt_min < pt < pt_max and 0 <= ipt < npt:
            if 2.0 < y < 2.5:
                deltaY = 2.5 - 2.0
                hist['sig_ypt1'][ipt] += SigTot * physicalWgt / (dpt * deltaY)
            if 2.5 < y < 3.0:
                deltaY = 3.0 - 2.5
                hist['sig_ypt2'][ipt] += SigTot * physicalWgt / (dpt * deltaY)
            if 3.0 < y < 3.5:
                deltaY = 3.5 - 3.0
                hist['sig_ypt3'][ipt] += SigTot * physicalWgt / (dpt * deltaY)
            if 3.5 < y < 4.0:
                deltaY = 4.0 - 3.5
                hist['sig_ypt4'][ipt] += SigTot * physicalWgt / (dpt * deltaY)
            if 4.0 < y < 4.5:
                deltaY = 4.5 - 4.0
                hist['sig_ypt5'][ipt] += SigTot * physicalWgt / (dpt * deltaY)
        if m_min < M < m_max and 0 <= im < nm:
            hist['sig_m'][im] += SigTot * physicalWgt / dm
            
            
    jac = (yVar_max - yVar_min)*(m2Var_max - m2Var_min)*(pt2Var_max - pt2Var_min)     
    
    vegasIntegrand = jac*SigTot

    return vegasIntegrand



