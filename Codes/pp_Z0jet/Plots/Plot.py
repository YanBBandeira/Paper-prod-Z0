import numpy as np
import matplotlib.pyplot as plt
import os



path = r'C:\Users\Callidus\Documents\Clones\Paper-prod-Z0\Codes\pp_Z0jet\Output\dsig_dpt.dat'
print("Existe?", os.path.exists(path))
print("Caminho absoluto:", os.path.abspath(path))

# ======================================================
# 1. Configurações gerais de estilo (equivalentes ao Gnuplot)
# ======================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "axes.labelsize": 22,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 18,
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "lines.linewidth": 3
})

# ======================================================
# 2. Caminhos dos arquivos (ajuste conforme necessário)
# ======================================================
output_dir = r'C:\Users\Callidus\Documents\Clones\Paper-prod-Z0\Codes\pp_Z0jet\Output'
data_dir =r'C:\Users\Callidus\Documents\Clones\Paper-prod-Z0\Codes\pp_Z0jet\DATA'

dist_Y = os.path.join(output_dir, "dsig_dy.dat")
dist_pT = os.path.join(output_dir, "dsig_dpt.dat")

ypT_files = {
    "2p25": os.path.join(output_dir, "dsig_dydpt_y2p25.dat"),
    "2p75": os.path.join(output_dir, "dsig_dydpt_y2p75.dat"),
    "3p25": os.path.join(output_dir, "dsig_dydpt_y3p25.dat"),
    "3p75": os.path.join(output_dir, "dsig_dydpt_y3p75.dat"),
    "4p25": os.path.join(output_dir, "dsig_dydpt_y4p25.dat")
}

data_files = {
    "Y": os.path.join(data_dir, "Ydist_dataset.csv"),
    "pT": os.path.join(data_dir, "pTdist_dataset.csv"),
    "2p25": os.path.join(data_dir, "YpTdist_2p25_dataset.csv"),
    "2p75": os.path.join(data_dir, "YpTdist_2p75_dataset.csv"),
    "3p25": os.path.join(data_dir, "YpTdist_3p25_dataset.csv"),
    "3p75": os.path.join(data_dir, "YpTdist_3p75_dataset.csv"),
    "4p25": os.path.join(data_dir, "YpTdist_4p25_dataset.csv")
}

save_dir = "plots"
os.makedirs(save_dir, exist_ok=True)

# ======================================================
# 3. Funções auxiliares
# ======================================================
def load_data(filename):
    """Carrega dados simples (duas colunas)"""
    return np.loadtxt(filename)

def load_dataset_with_errors(filename):
    """Carrega dados com erro (colunas: x, dx_low, y, dy_high)"""
    data = np.loadtxt(filename, skiprows=1)
    x, dx_low, y, dy_high = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    return x, y, dx_low, dy_high

# ======================================================
# 4. Função de plot genérica
# ======================================================
def make_plot(xlabel, ylabel, model_data, exp_data=None, title=None,
              logx=False, logy=False, label=None, xrange=None, yrange=None, filename=None):

    fig, ax = plt.subplots(figsize=(8, 6))

    # Escalas
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    # Model line
    ax.plot(model_data[:, 0], model_data[:, 1],
            color="blue", lw=3, label="IP-SAT")

    # Dados experimentais
    if exp_data is not None:
        x, y, ex, ey = exp_data
        ax.errorbar(x, y, xerr=np.abs(ex), yerr=np.abs(ey),
                    fmt="o", color="black", markersize=6,
                    label=r"LHCb data, $\sqrt{s} = 13$ TeV")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=20)
    if label:
        ax.text(0.85, 0.75, label, transform=ax.transAxes,
                fontsize=18, fontfamily="Helvetica")

    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.4)
    if xrange: ax.set_xlim(xrange)
    if yrange: ax.set_ylim(yrange)

    if filename:
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        print(f"✅ Saved: {filename}")
        plt.show()
    plt.close(fig)

# ======================================================
# 5. Plots individuais
# ======================================================

# --- pT distribution ---
model_pT = load_data(dist_pT)
exp_pT = load_dataset_with_errors(data_files["pT"])
make_plot(
    xlabel=r"$p_T$ [GeV]",
    ylabel=r"$d\sigma/dp_T$ [pb/GeV]",
    model_data=model_pT,
    exp_data=exp_pT,
    logx=True, logy=True,
    xrange=(1, 200),
    filename="pTdist.png"
)

# --- Y distribution ---
model_Y = load_data(dist_Y)
exp_Y = load_dataset_with_errors(data_files["Y"])
make_plot(
    xlabel=r"$Y$",
    ylabel=r"$d\sigma/dY$ [pb]",
    model_data=model_Y,
    exp_data=exp_Y,
    xrange=(2, 4.5),
    yrange=(0, 350),
    filename="Ydist.png"
)

# --- Y-pT distributions ---
y_labels = {
    "2p25": "2.0 < Y < 2.5",
    "2p75": "2.5 < Y < 3.0",
    "3p25": "3.0 < Y < 3.5",
    "3p75": "3.5 < Y < 4.0",
    "4p25": "4.0 < Y < 4.5"
}

for key, label in y_labels.items():
    model = load_data(ypT_files[key])
    exp = load_dataset_with_errors(data_files[key])
    make_plot(
        xlabel=r"$p_T$ [GeV]",
        ylabel=r"$d\sigma/dYdp_T$ [pb/GeV]",
        model_data=model,
        exp_data=exp,
        label=label,
        filename=f"YpTdist_{key}.png"
    )

print("\n✅ All plots saved in folder:", save_dir)
