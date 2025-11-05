#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include <Eigen/Dense>
#include "interpolation.h"  // cabeçalhos ALGLIB
#include <stdexcept>

namespace fs = std::filesystem;

// ======================================================
// 1. Physical constants and decay model
// ======================================================
class PhysicsParameters {
public:
    double pi;
    double alfem;
    double sin2;
    double aw;
    double Mz;
    double rs;
    double ml;

    PhysicsParameters() {
        pi = M_PI;
        alfem = 1.0 / 137.0;
        sin2 = 0.23;
        aw = std::asin(std::sqrt(sin2));
        Mz = 91.2;   // GeV
        rs = 13000.0; // sqrt(s)
        ml = 0.1056; // muon mass [GeV]
    }

    double dilepton_decay(double M) {
        double Mz2 = Mz * Mz;
        double M2 = M * M;
        double width = (alfem * M / (6.0 * std::pow(std::sin(2*aw), 2))) *
                       ((160.0/3.0)*std::pow(std::sin(aw), 4) - 40.0*std::pow(std::sin(aw), 2) + 21.0);
        double branch = 3.3 / 100.0;
        double inv_mass_dist = (1.0/pi) * ((M*width)/((M2 - Mz2)*(M2 - Mz2) + (M*width)*(M*width)));
        return inv_mass_dist * branch;
    }
};

// ======================================================
// 2. Histogram manager
// ======================================================
class HistogramManager {
public:
    std::vector<double> y_edges, pt_edges, m_edges;
    std::vector<double> sig_y, sig_pt, sig_m;
    std::vector<std::pair<double,double>> y_slices;
    std::vector<std::vector<double>> sig_ypt;

    HistogramManager(std::tuple<double,double,int> y_bins,
                     std::tuple<double,double,int> pt_bins,
                     std::tuple<double,double,int> m_bins)
    {
        auto [y_min, y_max, y_n] = y_bins;
        auto [pt_min, pt_max, pt_n] = pt_bins;
        auto [m_min, m_max, m_n] = m_bins;

        y_edges.resize(y_n);
        pt_edges.resize(pt_n);
        m_edges.resize(m_n);

        for (int i=0;i<y_n;i++) y_edges[i] = y_min + i*(y_max-y_min)/(y_n-1);
        for (int i=0;i<pt_n;i++) pt_edges[i] = pt_min + i*(pt_max-pt_min)/(pt_n-1);
        for (int i=0;i<m_n;i++) m_edges[i] = m_min + i*(m_max-m_min)/(m_n-1);

        sig_y.resize(y_n-1,0.0);
        sig_pt.resize(pt_n-1,0.0);
        sig_m.resize(m_n-1,0.0);

        y_slices = {{2.0,2.5},{2.5,3.0},{3.0,3.5},{3.5,4.0},{4.0,4.5}};
        sig_ypt.resize(y_slices.size(), std::vector<double>(pt_n-1,0.0));
    }

    void fill(double y, double pt, double M, double weight){
        // y distribution
        int iy = std::upper_bound(y_edges.begin(), y_edges.end(), y) - y_edges.begin() - 1;
        if (iy>=0 && iy<sig_y.size()) sig_y[iy]+=weight;

        // pt distribution
        int ipt = std::upper_bound(pt_edges.begin(), pt_edges.end(), pt) - pt_edges.begin() - 1;
        if (ipt>=0 && ipt<sig_pt.size()) sig_pt[ipt]+=weight;

        // M distribution
        int im = std::upper_bound(m_edges.begin(), m_edges.end(), M) - m_edges.begin() - 1;
        if (im>=0 && im<sig_m.size()) sig_m[im]+=weight;

        // y-pt double differential
        for (size_t idx=0; idx<y_slices.size(); idx++){
            double ymin=y_slices[idx].first, ymax=y_slices[idx].second;
            if (y>ymin && y<ymax && ipt>=0 && ipt<pt_edges.size()-1){
                double deltaY = ymax-ymin;
                sig_ypt[idx][ipt] += weight/deltaY;
            }
        }
    }

    void write_results(const std::string &outdir="Output"){
        fs::create_directory(outdir);

        auto write_vec = [&](const std::string &filename, const std::vector<double>& edges, const std::vector<double>& values){
            std::ofstream f(outdir+"/"+filename);
            for (size_t i=0;i<values.size();i++){
                double mid = 0.5*(edges[i+1]+edges[i]);
                f << mid << " " << values[i] << "\n";
            }
            f.close();
        };

        write_vec("dsig_dy.dat", y_edges, sig_y);
        write_vec("dsig_dpt.dat", pt_edges, sig_pt);
        write_vec("dsig_dm.dat", m_edges, sig_m);

        for (size_t i=0;i<y_slices.size();i++){
            std::ofstream f(outdir+"/dsig_dydpt_y"+std::to_string(y_slices[i].first)+"-"+std::to_string(y_slices[i].second)+".dat");
            for (size_t j=0;j<sig_ypt[i].size();j++){
                double mid = 0.5*(pt_edges[j+1]+pt_edges[j]);
                f << mid << " " << sig_ypt[i][j] << "\n";
            }
            f.close();
        }
        std::cout << "Histograms written to " << outdir << "/" << std::endl;
    }
};

// ======================================================
// 3. Grid interpolator
// ======================================================
class GridInterpolator {
public:
    int n_points;
    std::string filename;
    std::vector<double> y_grid, pt_grid, m_grid;
    std::vector<double> parton_values; // flattening 3D grid to 1D for ALGLIB
    alglib::spline2dinterpolant spline_y_pt;  // exemplo para 2D — veremos 3D abaixo

    GridInterpolator(const std::string &fname, int n=15)
      : filename(fname), n_points(n),
        y_grid(n_points), pt_grid(n_points), m_grid(n_points),
        parton_values(n_points*n_points*n_points, 0.0)
    {
        read_grid();
        build_interpolator();
    }

    void read_grid(){
        std::ifstream fin(filename);
        if(!fin.is_open()) throw std::runtime_error("Cannot open grid file");

        for(int i=0; i<n_points; i++){
            for(int j=0; j<n_points; j++){
                for(int k=0; k<n_points; k++){
                    double yv, ptv, mv, val;
                    fin >> yv >> ptv >> mv >> val;
                    y_grid[i] = yv;
                    pt_grid[j] = ptv;
                    m_grid[k] = mv;
                    parton_values[(i*n_points + j)*n_points + k] = val;
                }
            }
        }
        fin.close();
    }

    void build_interpolator(){
        //
        // Aqui está o *esqueleto* da lógica: ALGLIB não possui
        // diretamente “spline3dbuild” para grades retangulares regulares
        // documentadas (ao menos não no user‑guide público).
        // Uma abordagem prática: para cada fixed m_grid[k], construir
        // uma spline2d(y,pt) para os valores parton_grid[:,:,k], e depois
        // para um dado m, interpolar nas k vizinhas, e fazer interpolação linear em m.
        //

        // Exemplo simplificado para o nível m = m_grid[0] — repetir para todos k
        // e armazenar.
        //
        // Aqui, por simplicidade, faço *linha‐mestre* para m_fixed = m_grid[0]:
        std::vector<double> zy(n_points), zpt(n_points), zval(n_points);
        for(int i=0;i<n_points;i++){
            zy[i] = y_grid[i];
        }
        for(int j=0;j<n_points;j++){
            zpt[j] = pt_grid[j];
        }

        // Uso ALGLIB para cada k‑slice:
        for(int k=0; k<n_points; k++){
            // extrair valores f(i,j,k)
            std::vector<double> f2d(n_points * n_points);
            for(int i=0;i<n_points;i++){
                for(int j=0;j<n_points;j++){
                    f2d[i*n_points + j] = parton_values[(i*n_points + j)*n_points + k];
                }
            }
            alglib::spline2dbuildbilinearv(zy.data(), n_points,
                                          zpt.data(), n_points,
                                          f2d.data(), n_points,
                                          f2d.data(), n_points,
                                          spline2d[k]);
        }

        // Depois construir spline1d em m para cada (y,pt) célula, ou
        // usar interpolação linear entre slices.
    }

    double interpolate(double y, double pt, double m) const {
        //
        // Localizar k_low, k_high tal que m_grid[k_low] <= m <= m_grid[k_high]
        //
        int k_low = 0;
        while(k_low < n_points-1 && m > m_grid[k_low+1]) k_low++;
        if(k_low == n_points-1) k_low = n_points-2;
        int k_high = k_low + 1;

        double w = (m - m_grid[k_low]) / (m_grid[k_high] - m_grid[k_low]);

        // Avaliar spline2d para slice k_low e k_high
        double f_low  = alglib::spline2dcalc(spline2d[k_low], y, pt);
        double f_high = alglib::spline2dcalc(spline2d[k_high], y, pt);

        // Interpolação linear em m
        return (1.0 - w)*f_low + w*f_high;
    }

private:
    std::vector<alglib::spline2dinterpolant> spline2d;
};

// ======================================================
// 4. Cross section integrand (Fortran-style)
// ======================================================
class CrossSectionIntegrand {
public:
    PhysicsParameters *params;
    HistogramManager *hist;
    GridInterpolator *hadronic;

    CrossSectionIntegrand(PhysicsParameters *p, HistogramManager *h, GridInterpolator *had)
        : params(p), hist(h), hadronic(had) {}

    double operator()(const std::vector<double>& x){
        double yp = x[0], ym = x[1], ktp = x[2], ktm = x[3], phip = x[4], phim = x[5];
        double mp = params->ml;
        double mm = params->ml;

        // Kinematic cuts
        if(!((2.0<=yp && yp<=4.5) || (2.0<=ym && ym<=4.5))) return 0.0;
        if(!(ktp>=20.0 || ktm>=20.0)) return 0.0;

        // Transverse components
        double ktpx = ktp*std::cos(phip);
        double ktpy = ktp*std::sin(phip);
        double ktmx = ktm*std::cos(phim);
        double ktmy = ktm*std::sin(phim);

        double mperp_p = std::sqrt(ktp*ktp + mp*mp);
        double mperp_m = std::sqrt(ktm*ktm + mm*mm);

        double ptx = ktpx + ktmx;
        double pty = ktpy - ktmy; // Fortran-style
        double pt2 = ptx*ptx + pty*pty;
        double pt = std::sqrt(pt2);

        double deltaY = yp - ym;
        double m2 = mperp_p*mperp_p + mperp_m*mperp_m + 2.0*mperp_p*mperp_m*std::cosh(deltaY) - pt2;
        if(m2<=0) return 0.0;
        double M = std::sqrt(m2);

        double rs = params->rs;
        double xp = (ktp/rs)*std::exp(yp);
        double xm = (ktm/rs)*std::exp(ym);
        double xf = xp + xm;
        double y = std::log(xf * (rs / std::sqrt(pt2 + M*M)));

        double x1 = std::sqrt(M*M + pt2)/rs * std::exp(+y);
        double x2 = std::sqrt(M*M + pt2)/rs * std::exp(-y);

        if(!(2.0<=y && y<=4.5)) return 0.0;
        if(!(60.0<=M && M<=120.0)) return 0.0;
        if(x1>=1.0 || x2>=1.0) return 0.0;

        double varJacobian = (2.0/rs)*std::sqrt(M*M + pt2)*std::cosh(y);
        double preIntegral = (x1/(x1+x2))*varJacobian;

        double hadronic_val = hadronic->interpolate(y, pt, M);
        double decay = params->dilepton_decay(M);

        double sigma = preIntegral*decay*hadronic_val;

        hist->fill(y, pt, M, sigma);
        return sigma;
    }
};

// ======================================================
// 6. Main driver (VEGAS placeholder)
// ======================================================
int main(){
    PhysicsParameters params;
    HistogramManager hist({2.0,4.5,61}, {0.0,150.0,61}, {60.0,120.0,61});
    GridInterpolator hadronic("hadronic_grid.dat", 15);

    CrossSectionIntegrand integrand(&params, &hist, &hadronic);

    // Placeholder for VEGAS integration:
    std::cout << "VEGAS integration would run here (C++ version requires external library)\n";

    hist.write_results();
    return 0;
}
