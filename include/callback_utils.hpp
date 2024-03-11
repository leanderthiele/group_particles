/*! @file callback_utils.hpp
 *
 * @brief Collects some pre-implemented functionality that makes
 *        it easier to construct a subclass of the #Callback base class.
 *
 * By inheriting from some of the classes defined in the included
 * header files, the user can achieve a lot of functionality without
 * explicitly overriding the methods in #Callback.
 *
 * All classes are inside the #CallbackUtils namespace and take as their
 * first template argument a type constructed from the #AllFields template
 * (named AFields) which will not be explicitly documented for each of them.
 */

#ifndef CALLBACK_UTILS_HPP
#define CALLBACK_UTILS_HPP

#include "callback_utils_chunk.hpp"
#include "callback_utils_name.hpp"
#include "callback_utils_meta.hpp"
#include "callback_utils_meta_init.hpp"
#include "callback_utils_select.hpp"
#include "callback_utils_radius.hpp"
#include "callback_utils_grp_action.hpp"
#include "callback_utils_prt_action.hpp"
#include "callback_utils_prt_modify.hpp"

/*! @brief contains classes that implement parts of the #Callback base.
 *
 * By inheriting from classes in this namespace, the user can obtain a lot of functionality
 * that is fairly common quite easily, without the need to explicitly override methods
 * in the #Callback base.
 */
namespace CallbackUtils
{ }

#endif // CALLBACK_UTILS_HPP
