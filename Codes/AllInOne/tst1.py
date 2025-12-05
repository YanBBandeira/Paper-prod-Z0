import numpy as np
from scipy.integrate import dblquad # Substitui DADMUL/VEGAS para integração 2D
from scipy.interpolate import RegularGridInterpolator # Para FuncPartonLevelSigma

# --- 1. MÓDULOS DE CONSTANTES E GLOBAIS (Substituindo MODULEs) ---

class Parameters:
    PI = 4.0 * np.arctan(1.0)
    PI2 = PI**2.0
    ALFEM = 1.0 / 137.0
    SIN2 = 0.230  # Usando 0.23d0 do Fortran
    AW = np.arcsin(np.sqrt(SIN2))
    MZ = 91.2  # Massa do Z0 em GeV
    RS = 13000.0  # Energia CM sqrt(s) em GeV

class Globals:
    # Variáveis que seriam globais/common
    pt = 0.0  # ptZ fixo durante a integração interna
    x1 = 0.0
    x2 = 0.0
    m = 0.0

# Inicializa as classes de constantes e globais
params = Parameters()
globals_vars = Globals()

# --- 2. FUNÇÕES FÍSICAS (Simplificadas/Placeholders) ---

def DileptonDecay(M_var):
    """
    Calcula o termo de decaimento Z -> l+ l- (Baseado na função Fortran)
    """
    M = M_var
    M2 = M**2.0
    Mz2 = params.MZ**2.0

    # Termo de largura de decaimento (simplificado)
    DecayWidth = (params.ALFEM * M / (6.0 * (np.sin(2.0 * params.AW)**2.0))) * (
        (160.0 / 3.0) * (np.sin(params.AW)**4.0) - 40.0 * (np.sin(params.AW)**2.0) + 21.0
    )

    Branch = 3.3e-2  # 3.3%
    
    # Distribuição de massa invariante (Lorentzian)
    InvariantMassDist = (1.0 / params.PI) * (
        (M * DecayWidth) / ((M2 - Mz2)**2.0 + (M * DecayWidth)**2.0)
    )

    Result = InvariantMassDist * Branch
    return Result

# --- 3. SIMULAÇÃO DO GRID DE INTERPOLAÇÃO (FuncPartonLevelSigma) ---

# Como você mencionou que FuncPartonLevelSigma é um grid 2D (y, pt)
# Precisamos criar um grid de exemplo para demonstrar a interpolação.
# Em um código real, você carregaria estes dados de um arquivo.

# Exemplo de definição dos eixos do grid (deve corresponder aos seus bins)
NY_GRID = 60
NPT_GRID = 60
Y_MIN_GRID, Y_MAX_GRID = 2.0, 4.5
PT_MIN_GRID, PT_MAX_GRID = 0.0, 150.0

data = np.loadtxt(r"C:\Users\Callidus\Documents\Clones\Paper-prod-Z0\Codes\pp_Z0jet\Grids\DatFiles\tst_grid.dat")
y_axis = np.unique(data[:,0])
pt_axis = np.unique(data[:,1])
sigma_grid_data = data[:,2].reshape(len(y_axis ), len(pt_axis))
sigma_interpolator = RegularGridInterpolator((y_axis, pt_axis), sigma_grid_data, 
                                             bounds_error=False, fill_value=0.0)

def FuncPartonLevelSigma(yVar, ptVar, mVar = Parameters.MZ ):
    """
    Substitui a função Fortran pela interpolação do grid 2D (y, pt).
    O mVar (massa) é ignorado na interpolação 2D, mas mantido para consistência.
    """
    # O grid é 2D (y, pt). O mVar não é usado na interpolação 2D.
    sqrt_M2pT2 =  np.sqrt(mVar**2.0 + ptVar**2) 
    x1 = (sqrt_M2pT2/Parameters.RS)*np.exp(yVar)
    x2 = (sqrt_M2pT2/Parameters.RS)*np.exp(-yVar)
    VarJacobian = (2.0/Parameters.RS)* sqrt_M2pT2*np.cosh(yVar)
    preIntegral = x1/(x1 + x2)
    
    # Garantir que ptVar não seja negativo (embora o loop principal deva evitar isso)
    pt_interp = max(ptVar, PT_MIN_GRID) 
    
    # O interpolador espera um array de pontos: [[y1, pt1], [y2, pt2], ...]
    points = np.array([[yVar, pt_interp]])
    
    # O resultado é um array, pegamos o primeiro elemento [0]
    HadronicCrossSection = sigma_interpolator(points)[0] 

    # Em um código real, você calcularia x1, x2 aqui e usaria PDFs/TMDs
    # Aqui, retornamos o valor do grid.
    return HadronicCrossSection*VarJacobian*preIntegral

# --- 4. FUNÇÃO INTEGRANDA (Substituindo TestIntegrand) ---

def TestIntegrand(MZZ, yZ, ptZ):
    """
    Função integranda para DADMUL/dblquad.
    Integra sobre M e y, mantendo pT (ptZ) fixo.
    """
    # Variáveis de integração: MZZ (Massa), yZ (Rapidez)
    
    # Termo de decaimento
    DecayTerm = DileptonDecay(MZZ)
    
    # Seção de choque hadrônica (interpolação no grid)
    HadronicCS = FuncPartonLevelSigma(yZ, ptZ, MZZ)
    
    # O integrando final no Fortran era: HadronicCrossSection * 2.0 * Mzz
    # O fator 2*Mzz vem da transformação dM^2 -> 2M dM, se M^2 fosse a variável de integração.
    # Como estamos integrando sobre M (MZZ), o fator 2*Mzz é o jacobiano de M^2 -> M.
    
    Result = DecayTerm * HadronicCS 
    
    return Result



# --- 5. PROGRAMA PRINCIPAL (Substituindo z_pt_distribution) ---

def z_pt_distribution_python():
    

    # --- Binning para pT ---
    PT_MIN = np.log10(1.0)
    PT_MAX = np.log10(400.0)
    N_POINTS = 80
    DPT = (PT_MAX - PT_MIN) / N_POINTS
    
    GEV_TO_PB = 0.389e9 # GeV^-2 to pb
    
    # --- Loop em pT ---
    
    ptZ_values = []
    dsigma_values = []
    
    print("\n--- Starting pT Loop ---")
    
    for ipt in range(1, N_POINTS + 1):
        # Cálculo do ptZ (logarítmico)
        ptZ_log = PT_MIN + (ipt - 1.0) * DPT
        ptZ = 10.0**ptZ_log  # Ajuste de offset como no Fortran
        
        # Fixa a variável global pt para ser usada no integrando
        globals_vars.pt = ptZ 
        
        # --- Integração 2D (Substituindo DADMUL) ---
        # dblquad(func, a, b, gfun, hfun) -> func(y, x)
        # Aqui, func(MZZ, yZ) onde ptZ é fixo.
        
        # Limites de integração para M e y (do Fortran)
        M_MIN, M_MAX = 60.0**2.0, 120.0**2.0
        Y_MIN, Y_MAX = 2.0, 4.5
        
        # A função TestIntegrand espera (MZZ, yZ) na ordem que dblquad chama.
        # dblquad(func, a, b, gfun, hfun) -> func(x, y) onde x é a variável externa (MZZ)
        # e y é a variável interna (yZ).
        
        # Definindo os limites internos (yZ) como função da variável externa (MZZ)
        # Como os limites de y e M são constantes, gfun e hfun são constantes.
        
        # Definindo a função lambda para chamar TestIntegrand na ordem correta (MZZ, yZ)
        # dblquad chama func(x, y) onde x é a variável externa (M_MAX/MIN) e y a interna (Y_MAX/MIN)
        
        # Vamos definir a ordem de integração: M (externa), y (interna)
        # func(MZZ, yZ)
        
        # Definindo a função lambda para garantir a ordem de chamada correta para dblquad
        # dblquad(func, a, b, gfun, hfun) -> func(x, y)
        # x (externo) = MZZ, y (interno) = yZ
        
        integrand_wrapper = lambda yZ, MZZ: TestIntegrand(MZZ, yZ, ptZ)
        
        # Executa a integração 2D
        # Resultado: (Integral, Erro Estimado)
        Result_integral, error = dblquad(
            integrand_wrapper, 
            M_MIN, M_MAX,  # Limites externos (MZZ)
            lambda MZZ: Y_MIN, lambda MZZ: Y_MAX, # Limites internos (yZ)
            epsabs=1e-5, epsrel=1e-5 # Precisão similar ao EPS do Fortran
        )
        
        dsigma = GEV_TO_PB * Result_integral
        
        ptZ_values.append(ptZ)
        dsigma_values.append(dsigma)
        
        print(f"pT={ptZ:.4f} GeV, dsigma/dpT = {dsigma:.4e} pb")

    # --- Saída de Dados ---
    print("\n--- Output ---")
    with open('dsig_dpt.dat', 'w') as f:
        for pt, dsigma in zip(ptZ_values, dsigma_values):
            f.write(f"{pt:10.4f}  {dsigma:12.5e}\n")
    
    print("Calculation complete. Results saved to dsig_dpt.dat")





# Executar o programa principal
if __name__ == "__main__":
    z_pt_distribution_python()
    
    
import numpy as np
import matplotlib.pyplot as plt

# Carregar dados
data = np.loadtxt("dsig_dpt.dat")

pt = data[:,0]
dsig = data[:,1]



pt_low  = np.array([0.0, 2.2, 3.4, 4.6, 5.8, 7.2, 8.7, 10.5, 12.8, 15.4, 19.0, 24.5, 34.0, 63.0])
pt_high = np.array([2.2, 3.4, 4.6, 5.8, 7.2, 8.7, 10.5, 12.8, 15.4, 19.0, 24.5, 34.0, 63.0, 270.0])

pt_bins = 0.5*(pt_low + pt_high)  # centro do bin
sigma_pt = np.array([5.70, 11.07, 11.44, 11.28, 9.94, 8.86, 7.75, 6.44,
                     5.16, 4.03, 2.88, 1.774, 0.674, 0.0361])
err_pt = np.array([0.5, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2])



# Plot
plt.figure(figsize=(8,6))
plt.plot(pt, dsig, marker='', linestyle='-', color='navy', label=r"$\frac{d\sigma}{dp_T}$")
plt.plot(pt_bins, sigma_pt, 'o', color='red', label='LHCb')
plt.xlabel(r"$p_T$  [GeV]", fontsize=14)
plt.ylabel(r"$\frac{d\sigma}{dp_T}$  [pb/GeV]", fontsize=14)
plt.title("Distribuição de $p_T$ do $Z^0$", fontsize=15)

#plt.yscale('log')
plt.xscale("log")      # Geralmente faz sentido para distribuição de pT
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend(fontsize=13)

plt.tight_layout()
plt.show()



