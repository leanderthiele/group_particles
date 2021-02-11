#ifndef HALO_PARTICLES_HPP
#define HALO_PARTICLES_HPP

#include <functional>
#include <cstdio>

#include "fields.hpp"
#include "illustris_fields.hpp"
#include "callback.hpp"
#include "callback_utils.hpp"
#include "read_fields.hpp"
#include "workspace.hpp"
#include "workspace_memory.hpp"
#include "grp_loop.hpp"
#include "prt_loop.hpp"

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
    template_checks<GroupFields, ParticleFields>();

    Workspace<GroupFields,ParticleFields> ws;

    ws.grp_loop(callback);

    ws.prt_loop(callback);
}

#endif // HALO_PARTICLES_HPP
