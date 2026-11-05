// P2-2d (a3)/(a5): dump BOTH toughness evaluators from the REAL code.
//
// ⚠⚠ THE WHOLE POINT IS THAT THESE ARE THE CODE'S OWN EVALUATORS. Rev 2 of the
// packet read 34 % closure from a HAND interpolation taken at two different
// temperatures. `(a3)` therefore requires `KIcNasseri` and
// `Grimsel::K_Ic_of` themselves, compiled from src/, never a transcription.
#include <cstdio>
#include "Numeric/SpCriterion.H"
#include "Numeric/Material/Grimsel.H"
int main()
{
    Numeric::Material::Grimsel g;
    printf("T_K THETA_C KIC_NASSERI KIC_GRIMSEL E_PA STRAIN\n");
    for (double T = 293.15; T <= 1223.15 + 1e-9; T += 1.0)
        printf("%.4f %.4f %.17e %.17e %.17e %.17e\n",
               T, T - 273.15,
               Numeric::SpCriterion::KIcNasseri(T, 1.0),
               g.K_Ic_of(T), g.E_of(T), g.strain_of(T));
    return 0;
}
