#include "rhoThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "rhoConst.H"

#include "hConstThermo.H"
//#include "penaltyHConstThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "porousHeRhoThermo.H"
#include "pureMixture.H"


namespace Foam{
makeThermos
(
    rhoThermo,
    porousHeRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    rhoConst,
    specie
);

makeThermos
(
    rhoThermo,
    porousHeRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);
}
