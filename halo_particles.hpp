#ifndef HALO_PARTICLES_HPP
#define HALO_PARTICLES_HPP

#include <functional>
#include <cstdio>

#include "fields.hpp"
#include "illustris_fields.hpp"

// Main function to be used in the interface
// 
//
// halo_catalog and particle_catalog should be regular expressions telling the
// function where to find the simulation data.
//
// 
template<uint8_t PartType, typename GroupFields, typename ParticleFields>
void
halo_particles (const std::string &halo_catalog,
                const std::string &particle_catalog,
                std::function<bool(void **)> halo_selector,
                std::function<float(void **)> halo_radius,
                std::function<bool(void **, void **)> particle_selector,
                std::function<void(void **, void **)> particle_action)
{
    template_checks<GrpFields, ParticleFields>();
}

#endif // HALO_PARTICLES_HPP
