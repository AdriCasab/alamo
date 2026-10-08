#include <cstdio>
#include <cmath>
#include "g_standalone.H"
int main(){
    Numeric::Material::Grimsel g;
    // dump the REAL C++ strain at the temperatures (b6) spans
    printf("THETA_C STRAIN E_PA SIGMA_PA\n");
    for (double T=293.15; T<=1000.15; T+=1.0) {
        printf("%.4f %.17e %.17e %.17e\n",
               T-273.15, g.strain_of(T), g.E_of(T), g.sigma_thermal(T,293.15));
    }
    return 0;
}
