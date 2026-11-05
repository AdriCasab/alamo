// P2-2b (b2): dump the REAL C++ law so an offline twin can score it.
//
// ⚠ the REAL law, straight from src/ -- never a copy. See standalone/README.md
#include <cstdio>
#include <cmath>
#include "Numeric/Material/Grimsel.H"

// ⚠⚠ TWO SWEEPS, AND THE SECOND IS THE ONE WITH POWER.
//
//   NODE  708 points at INTEGER theta, 20..727 degC -- the packet's (b2) sweep.
//   OFF   708 points at theta = 20 + i*(980/707), which is irrational to the
//         1 K grid and so lands BETWEEN nodes almost everywhere.
//
// P2-2b made E_of a dense uniform-grid lookup. Its correctness claim is that
// linear interpolation between 1 K nodes is EXACT because E/E_0 is piecewise
// linear with every Table 3 kink on an integer degree. A sweep evaluated only
// at integer theta hits only grid nodes, where ANY interpolation scheme agrees
// by construction -- it could not detect a kink-straddling error, which is the
// one way this change could be wrong. The OFF sweep is what tests it.
int main()
{
    Numeric::Material::Grimsel g;
    printf("SWEEP THETA_C STRAIN E_PA SIGMA_PA\n");
    for (int i = 0; i < 708; ++i)
    {
        const double th = 20.0 + double(i);
        const double T  = th + 273.15;
        printf("NODE %.17e %.17e %.17e %.17e\n",
               th, g.strain_of(T), g.E_of(T), g.sigma_thermal(T, 293.15));
    }
    for (int i = 0; i < 708; ++i)
    {
        const double th = 20.0 + double(i) * (980.0 / 707.0);
        const double T  = th + 273.15;
        printf("OFF %.17e %.17e %.17e %.17e\n",
               th, g.strain_of(T), g.E_of(T), g.sigma_thermal(T, 293.15));
    }
    return 0;
}
