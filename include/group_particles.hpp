/*! @file group_particles.hpp
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

/*! @mainpage Overview
 *
 * The code is header-only, which can be accessed by including the file group_particles.hpp
 * in the user's source code.
 *
 * The single exposed routine is #group_particles, with the signature
 * @code
 * template<typename AFields>
 * void group_particles (Callback<AFields> &callback);
 * @endcode
 * This routine will load the group and particle catalog(s) from disk and perform the user defined
 * actions on them.
 *
 * The template parameter `AFields` defines which data fields from the group and particle catalogs
 * should be loaded into memory and made accessible.
 * The `AFields` type should be constructed using the #AllFields template, with signature
 * @code
 * template<typename GroupFields, typename ParticleFields>
 * struct AllFields;
 * @endcode
 * Here, `GroupFields` and `ParticleFields` should be constructed using the `GrpFields` and `PrtFields`
 * types. These have signatures
 * @code
 * template<typename... Fields>
 * struct GrpFields;
 *
 * template<typename... Fields>
 * struct PrtFields;
 * @endcode
 * The template parameter packs are chosen by the user. The types passed can be constructed using the
 * #FIELD macro; a number of examples for Illustris- and Gadget-type simulations are provided in
 * common_fields.hpp.
 * 
 * Having specified which data fields are to be read from the simulation files, we now need to specify
 * which actions are to be performed on them.
 * To this end, the user should inherit from the #Callback class and override the methods defined there.
 * Thus, we obtain an implemented #Callback-subclass which can then be passed to the #group_particles
 * routine.
 * 
 * Some of the methods in the #Callback class are simply meant to inform the code of the layout of the
 * data files, namely
 *      - #Callback::grp_chunk, #Callback::prt_chunk for the file names,
 *      - #Callback::grp_name, #Callback::prt_name for the internal file layout,
 *      - #Callback::read_grp_meta, #Callback::read_prt_meta for some meta-data the code requires.
 *
 * The code also allows the user to initially read some meta-data (probably cosmological parameters,
 * mass tables, etc.) from the data files and store them in their own #Callback subclass instance.
 * These routines are
 *      - #Callback::read_grp_meta_init, #Callback::read_prt_meta_init.
 *
 * Now, we need to define specifically what to do wih the groups and particles.
 * This is accomplished by overriding the following methods:
 *      - #Callback::grp_select defines which groups should be considered,
 *      - #Callback::grp_action lets the user do something with a group's data
 *                              (typically store some group properties for later use),
 *      - #Callback::grp_radius defines how to compute a group's radius,
 *      - #Callback::prt_action defines what to do with the particles that fall into a group's radius.
 *
 * All these methods take at least one of the two auxiliary types
 * #Callback::GrpProperties and #Callback::PrtProperties.
 * These types contain the data fields defined by the `AFields` type above.
 * From the user's side, the only interface that should be used to these types is the `get` template method
 * (#Callback::BaseProperties::get).
 * It allows to retrieve individual properties of a group or particle.
 *
 * Because many applications will require very similar implementations of many of the #Callback methods,
 * we provide a number of them in the #CallbackUtils namespace.
 *
 */

/*! @brief Runs the code.
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
group_particles (Callback<AFields> &callback)
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
