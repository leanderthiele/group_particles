#ifndef HALO_PARTICLES_HPP
#define HALO_PARTICLES_HPP

#include <functional>
#include <cstdio>

#include "fields.hpp"
#include "illustris_fields.hpp"
#include "callback.hpp"

// Main function to be used in the interface
// 
//
// halo_catalog and particle_catalog should be regular expressions telling the
// function where to find the simulation data.
//
// 
template<typename GroupFields, typename ParticleFields>
void
halo_particles (Callback &callback)
{
    template_checks<GrpFields, ParticleFields>();
}

#endif // HALO_PARTICLES_HPP
