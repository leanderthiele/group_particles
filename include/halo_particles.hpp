/*! @file halo_particles.hpp
 *
 * @brief For convenience, the user can just include this single header file.
 *        It contains the single routine to be called.
 */

#ifndef HALO_PARTICLES_HPP
#define HALO_PARTICLES_HPP

#include "fields.hpp"
#include "common_fields.hpp"
#include "callback.hpp"
#include "callback_utils.hpp"
#include "hdf5_utils.hpp"
#include "workspace.hpp"
#include "workspace_memory.hpp"
#include "workspace_sorting.hpp"
#include "workspace_meta_init.hpp"
#include "grp_loop.hpp"
#include "prt_loop.hpp"

/*! Runs the code.
 *
 * @tparam AFields      a type constructed from the #AllFields template,
 *                      defines which fields the code should read from the data files.
 * @param[in,out] callback      by subclassing from the #Callback abstract base class,
 *                              the user defines which functionality they want the
 *                              code to fulfill.
 *                              Furthermore, the passed instance can be used to store
 *                              data.
 *
 * See the documentation of the #Callback class for all methods that need to be overriden.
 * Here, we give the order in which the non-const member functions are called
 * (so that the user knows when data that they want to store inside the passed instance
 * will be available).
 *
 * 1. allow the user to read meta-data from the first group and particle files:
 *    - #Callback::read_grp_meta_init
 *    - #Callback::read_prt_meta_init
 *
 * 2. perform user-defined action on each group:
 *    - #Callback::grp_action
 *
 *    Will be called consecutively for each group that passes #Callback::grp_select.
 *
 * 3. perform user-defined action on each particle in each group:
 *    - #Callback::prt_action
 *
 *    Will be called in an undefined order, perhaps in parallel.
 *    However, it is guaranteed that no two calls will be simultaneous if they pass
 *    the same group.
 *    Thus, if the user's #Callback::prt_action implementation is local
 *    (in the sense that it only acts on data associated with a single group;
 *     this will be the case in almost all applications)
 *    the user is not required to take any precautations with regard to thread safety.
 */
template<typename AFields>
void
halo_particles (Callback<AFields> &callback)
{
    #ifndef NDEBUG
    AFields::print_field_info();
    #endif // NDEBUG

    Workspace<AFields> ws { callback };

    ws.meta_init();

    ws.grp_loop();

    ws.prt_loop();
}

#endif // HALO_PARTICLES_HPP
