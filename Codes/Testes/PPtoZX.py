import numpy as np
import vegas
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

# -----------------------------------------------------
# Constantes físicas
# -----------------------------------------------------
pi = np.pi
alfem = 1.0 / 137.0
sin2 = 0.231
aw = np.arcsin(np.sqrt(sin2))
Mz = 91.2
rs = 13000.0
mp = mm = 0.1056

# -----------------------------------------------------
# Leitura do grid parton-level
# -----------------------------------------------------
data = np.loadtxt(r"C:\Users\Callidus\Documents\Clones\Paper-prod-Z0\Codes\Testes\tst_grid.dat")
y_grid = np.unique(data[:,0])
pt_grid = np.unique(data[:,1])
f_grid = data[:,2].reshape(len(y_grid), len(pt_grid))
interp = RegularGridInterpolator((y_grid, pt_grid), f_grid,
                                 bounds_error=False, fill_value=0.0)

# -----------------------------------------------------
# Função de decaimento dileptônico
# -----------------------------------------------------
def dilepton_decay(M):
    M2 = M*M
    Mz2 = Mz*Mz
    decay_width = (alfem * M) / (6.0 * (np.sin(2.0 * aw)**2)) * (
        (160.0/3.0)*np.sin(aw)**4 - 40.0*(np.sin(aw)**2) + 21.0
    )
    branch = 3.366 / 100.0
    inv_mass_dist = (1.0/pi) * (M * decay_width) / ((M2 - Mz2)**2 + (M * decay_width)**2)


    return inv_mass_dist * branch 

# -----------------------------------------------------
# Função integrando (VEGAS)
# -----------------------------------------------------
def integrand(x):
    # Cada componente já é float
    x0, x1, x2, x3, x4, x5 = x

    # Variáveis físicas
    yp = 2.0 + (4.5 - 2.0) * x0
    ym = 2.0 + (4.5 - 2.0) * x1
    ktp = 20.0 + (200.0 - 20.0) * x2
    ktm = 20.0 + (200.0 - 20.0) * x3
    phip = pi*x4
    phim = pi*x5

    # Momento dos múons
    pxp = ktp*np.cos(phip)
    pyp = ktp*np.sin(phip)
    pzp = ktp*np.sinh(yp)
    Ep  = np.sqrt(mp**2 + pxp**2 + pyp**2 + pzp**2)

    pxm = ktm*np.cos(phim)
    pym = ktm*np.sin(phim)
    pzm = ktm*np.sinh(ym)
    Em  = np.sqrt(mm**2 + pxm**2 + pym**2 + pzm**2)

    # Z boson
    pxZ = pxp + pxm
    pyZ = pyp + pym
    pzZ = pzp + pzm
    ptZ = np.sqrt(pxZ**2 + pyZ**2)
    EZ  = Ep + Em
    M2  = EZ**2 - (pxZ**2 + pyZ**2 + pzZ**2)
    if M2 <= 0:
        return 0.0
    M = np.sqrt(M2)
    yZ = 0.5*np.log((EZ+pzZ)/(EZ-pzZ))

    # Cortes de aceitação
    eta_p = 0.5*np.log((Ep+pzp)/(Ep-pzp))
    eta_m = 0.5*np.log((Em+pzm)/(Em-pzm))
    if ktp < 20 or ktm < 20 or eta_p < 2 or eta_p > 4.5 or eta_m < 2 or eta_m > 4.5:
        return 0.0

    # Jacobiano
    jac_phase = ((4.5-2.0)**2 * (200.0-20.0)**2 * (2*pi*ktp) * (2*pi*ktm))
    var_jacobian = (2/rs) * np.sqrt(M**2 + ptZ**2) * np.cosh(yZ)

    hadronic_cs = interp([yZ, ptZ])[0] * 2.568e-9 # pb -> GeV^-2
    GeV2_to_pb = 0.389379e9  # 1 GeV^-2 = 0.389379e9 pb
    return (dilepton_decay(M)  * hadronic_cs * jac_phase * var_jacobian * GeV2_to_pb) / M2

# -----------------------------------------------------
# Rodar VEGAS
# -----------------------------------------------------
integ = vegas.Integrator([[0,1]]*6)
result = integ(integrand, nitn=10, neval=20000)
print(result.summary())

# -----------------------------------------------------
# Gerar amostras ponderadas usando VEGAS
# -----------------------------------------------------
N = 5000
ys, pts, ms, weights = [], [], [], []

for x, wgt in integ.random():  # pega N pontos do integrador
    x0, x1, x2, x3, x4, x5 = x

    yp = 2.0 + (4.5 - 2.0) * x0
    ym = 2.0 + (4.5 - 2.0) * x1
    ktp = 20.0 + (200.0 - 20.0) * x2
    ktm = 20.0 + (200.0 - 20.0) * x3
    phip = pi*x4
    phim = pi*x5

    pxp = ktp*np.cos(phip)
    pyp = ktp*np.sin(phip)
    pzp = ktp*np.sinh(yp)
    Ep  = np.sqrt(mp**2 + pxp**2 + pyp**2 + pzp**2)

    pxm = ktm*np.cos(phim)
    pym = ktm*np.sin(phim)
    pzm = ktm*np.sinh(ym)
    Em  = np.sqrt(mm**2 + pxm**2 + pym**2 + pzm**2)

    pxZ = pxp + pxm
    pyZ = pyp + pym
    pzZ = pzp + pzm
    ptZ = np.sqrt(pxZ**2 + pyZ**2)
    EZ  = Ep + Em
    M2  = EZ**2 - (pxZ**2 + pyZ**2 + pzZ**2)
    if M2 <= 0:
        continue
    M = np.sqrt(M2)
    yZ = 0.5*np.log((EZ+pzZ)/(EZ-pzZ))

    eta_p = 0.5*np.log((Ep+pzp)/(Ep-pzp))
    eta_m = 0.5*np.log((Em+pzm)/(Em-pzm))
    if ktp < 20 or ktm < 20 or eta_p < 2 or eta_p > 4.5 or eta_m < 2 or eta_m > 4.5:
        continue

    jac_phase = ((4.5-2.0)**2 * (200.0-20.0)**2 * (2*pi*ktp) * (2*pi*ktm))
    var_jacobian = (2/rs) * np.sqrt(M**2 + ptZ**2) * np.cosh(yZ)
    hadronic_cs  = interp([yZ, ptZ])[0] * 2.568e-9 # pb -> GeV^-2
    GeV2_to_pb = 0.389379e9  # 1 GeV^-2 = 0.389379e9 pb

    w = (dilepton_decay(M)  * hadronic_cs * jac_phase * var_jacobian * GeV2_to_pb) / (M2)

    ys.append(yZ)
    pts.append(ptZ)
    ms.append(M)
    weights.append(w)

ys = np.array(ys)
pts = np.array(pts)
ms = np.array(ms)
weights = np.array(weights)







# Exemplo Tabela 14
y_low  = np.array([2.000, 2.125, 2.250, 2.375, 2.500, 2.625, 2.750, 2.875,
                   3.000, 3.125, 3.250, 3.375, 3.500, 3.625, 3.750, 3.875,
                   4.000, 4.250])
y_high = np.array([2.125, 2.250, 2.375, 2.500, 2.625, 2.750, 2.875, 3.000,
                   3.125, 3.250, 3.375, 3.500, 3.625, 3.750, 3.875, 4.000,
                   4.250, 4.500])
y_bins = 0.5*(y_low + y_high)  # centro do bin
sigma_y = np.array([12.8, 40.4, 65.2, 87.5, 106.3, 122.7, 134.5, 141.7,
                    147.5, 145.4, 134.8, 118.5, 99.0, 77.6, 57.9, 39.5,
                    18.2, 2.7])        # dσ/dy
err_y = np.array([1.0, 0.9, 0.9, 0.8, 0.8, 0.7, 0.7, 0.6, 0.6, 0.5])                     # erro total (exemplo)



# -----------------------------------------------------
# Histogramas
# -----------------------------------------------------
plt.figure()
# Dados LHCb
plt.plot(y_bins, sigma_y, 'o-', color='red', label='LHCb')
# Teoria (do VEGAS)
plt.hist(ys, bins=50, weights=weights, density=False, histtype='step', color='blue', label='Teoria')
plt.xlabel("y_Z")
plt.ylabel("dσ/dy [pb]")
plt.title("Figura 3: dσ/dy_Z")
plt.legend()
plt.show()


# Exemplo: extração da Tabela 15
pt_low  = np.array([0.0, 2.2, 3.4, 4.6, 5.8, 7.2, 8.7, 10.5, 12.8, 15.4, 19.0, 24.5, 34.0, 63.0])
pt_high = np.array([2.2, 3.4, 4.6, 5.8, 7.2, 8.7, 10.5, 12.8, 15.4, 19.0, 24.5, 34.0, 63.0, 270.0])

pt_bins = 0.5*(pt_low + pt_high)  # centro do bin
sigma_pt = np.array([5.70, 11.07, 11.44, 11.28, 9.94, 8.86, 7.75, 6.44,
                     5.16, 4.03, 2.88, 1.774, 0.674, 0.0361])
err_pt = np.array([0.5, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2])

plt.figure()
plt.plot(pt_bins, sigma_pt, 'o-', color='red', label='LHCb')
plt.hist(pts, bins=50, weights=weights, density=False, histtype='step', color='blue', label='Teoria')
plt.xlabel("pT_Z [GeV]")
plt.ylabel("dσ/dpT [pb/GeV]")
plt.title("Figura 4: dσ/dpT_Z")
plt.xscale("log")
plt.legend()
plt.show()


plt.figure()
plt.hist(ms, bins=40, weights=weights, density=True)
plt.xlabel("M [GeV]")
plt.ylabel("dσ/dM")
plt.show()

plt.figure()
plt.hist2d(pts, ys, bins=(40,40), weights=weights, density=True, cmap='viridis')
plt.xlabel("pT_Z [GeV]")
plt.ylabel("y_Z")
plt.title("d²σ/(dpT dy)")
plt.colorbar(label="arbitrary units")
plt.show()
