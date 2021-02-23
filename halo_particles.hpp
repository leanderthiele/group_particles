#ifndef HALO_PARTICLES_HPP
#define HALO_PARTICLES_HPP

#include <functional>
#include <cstdio>

#include "fields.hpp"
#include "illustris_fields.hpp"
#include "callback.hpp"
#include "callback_utils.hpp"
#include "hdf5_utils.hpp"
#include "workspace.hpp"
#include "workspace_memory.hpp"
#include "workspace_sorting.hpp"
#include "grp_loop.hpp"
#include "prt_loop.hpp"

// Main function to be used in the interface
template<typename AFields>
void
halo_particles (Callback<AFields> &callback)
{
    #ifndef NDEBUG
    AFields::print_field_info();
    #endif // NDEBUG

    Workspace<AFields> ws { callback };

    ws.grp_loop();

    ws.prt_loop();
}

#endif // HALO_PARTICLES_HPP
