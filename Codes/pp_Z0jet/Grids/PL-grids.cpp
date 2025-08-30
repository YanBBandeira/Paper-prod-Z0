#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

// Constantes globais
namespace parameters {
    constexpr double pi = 4.0 * atan(1.0);
    constexpr double pi2 = pi * pi;
    constexpr double alfem = 1.0 / 137.0;
    constexpr double sin2 = 0.23;
    constexpr double aw = asin(sqrt(sin2));
    constexpr double Mz = 91.2; // Z0
    constexpr double rs = 13000.0; // sqrt(s)
}

// Prototipos
double funcPartonLevelSigma(double ptVar, double yVar);
double IntegrandHadronicCrossSection(double alf);
void PartonTargetCrossSection(double &Fvar, double ptVar, double zVar, double mVar, double mfVar, double gfv, double gfa);
double TransverseMomentumIntegral(int N, const std::vector<double>& X);
double TransverseMomentumIntegrand(double akt, double atheta);
double Epsilon1(double akt, double atheta, double pt, double z, double epps);
double Epsilon2(double akt, double atheta, double pt, double z, double epps);
double ugd(double akt);
double F_KS(double x2, double akt);
double sc(double Q2);

// Função principal
int main() {
    constexpr int nPoints = 5000;
    std::vector<double> y(nPoints), pt(nPoints);
    double y_min = 1.0, y_max = 5.0;
    double pt_min = log10(1.0), pt_max = log10(150.0);
    double dy = (y_max - y_min) / (nPoints - 1);
    double dpt = (pt_max - pt_min) / (nPoints - 1);

    std::ofstream fout("kslinear_grid.dat");

    for (int iy = 0; iy < nPoints; ++iy) {
        y[iy] = y_min + iy * dy;
        for (int ipt = 0; ipt < nPoints; ++ipt) {
            pt[ipt] = pow(10.0, pt_min + ipt * dpt) - 0.9;
            double partonLevelSigma = funcPartonLevelSigma(pt[ipt], y[iy]);
            std::cout << "Computing point: y = " << y[iy] << " pt = " << pt[ipt] << std::endl;
            fout << y[iy] << " " << pt[ipt] << " " << partonLevelSigma << std::endl;
        }
    }
    fout.close();
    return 0;
}

// Implementações simplificadas (algumas dependem de rotinas externas)
double funcPartonLevelSigma(double ptVar, double yVar) {
    using namespace parameters;
    double M2 = Mz * Mz;
    double pt2 = ptVar * ptVar;
    double x1 = (sqrt(M2 + pt2) / rs) * exp(yVar);
    double x2 = (sqrt(M2 + pt2) / rs) * exp(-yVar);

    // Aqui você deve implementar a integração numérica (ex: GSL)
    double result = 0.0; // dgauss(IntegrandHadronicCrossSection, 0.01, 1.0, 1e-4);
    double units = 0.389e9; // GeV^-2 to pb
    return result * units;
}

// As demais funções devem ser convertidas de forma similar.
// Para integração multidimensional, use GSL ou outra biblioteca.
// Para PDFs, use LHAPDF ou TMDlib (há bindings em C++).

// Exemplo de stub para sc:
double sc(double Q2) {
    using namespace parameters;
    double lambqcd = 0.35;
    double lambqcd2 = lambqcd * lambqcd;
    double b0 = 11.0 - (10.0 / 3.0);
    return (4.0 * pi) / (b0 * log(Q2 / lambqcd2));
}

// As funções Epsilon1, Epsilon2, ugd, F_KS, etc. podem ser convertidas diretamente.
// As rotinas de integração e PDFs exigem bibliotecas externas.
