#include <cmath>
#include <iostream>

namespace parameters {
    constexpr double pi = 4.0 * atan(1.0);
    // other parameters can be added here if needed
}

// Prototype / stub for the physics routine (implement your logic there)
void CrossSectionDimuon(double phiP, double phiM,
                        double kP, double kM,
                        double yP, double yM,
                        double &dsigma) {
    // ...existing code...
    dsigma = 0.0; // replace with real calculation
}

int main() {
    using namespace parameters;

    const int nPoints = 100; // adjust as needed

    // ranges
    double phiPmax = pi, phiPmin = 0.0;
    double phiMmax = pi, phiMmin = 0.0;

    double kPmax = 200.0, kPmin = 20.0;
    double kMmax = 200.0, kMmin = 20.0;

    double yPmax = 4.5, yPmin = 2.0;
    double yMmax = 4.5, yMmin = 2.0;

    // step sizes (use nPoints as number of subdivisions)
    double dphip = (phiPmax - phiPmin) / nPoints;
    double dphim = (phiMmax - phiMmin) / nPoints;
    double dkP   = (kPmax   - kPmin)   / nPoints;
    double dkM   = (kMmax   - kMmin)   / nPoints;
    double dyP   = (yPmax   - yPmin)   / nPoints;
    double dyM   = (yMmax   - yMmin)   / nPoints;

    double dsigma = 0.0;

    for (int iphiP = 0; iphiP <= nPoints; ++iphiP) {
        double phiP = phiPmin + iphiP * dphip;
        for (int iphiM = 0; iphiM <= nPoints; ++iphiM) {
            double phiM = phiMmin + iphiM * dphim;
            for (int ikP = 0; ikP <= nPoints; ++ikP) {
                double kP = kPmin + ikP * dkP;
                for (int ikM = 0; ikM <= nPoints; ++ikM) {
                    double kM = kMmin + ikM * dkM;
                    for (int iyP = 0; iyP <= nPoints; ++iyP) {
                        double yP = yPmin + iyP * dyP;
                        for (int iyM = 0; iyM <= nPoints; ++iyM) {
                            double yM = yMmin + iyM * dyM;
                            CrossSectionDimuon(phiP, phiM, kP, kM, yP, yM, dsigma);
                            // use dsigma (accumulate, write out, etc.)
                        }
                    }
                }
            }
        }
    }

    return 0;
}

