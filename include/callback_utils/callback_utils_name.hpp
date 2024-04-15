/*! @file callback_utils_name.hpp
 *
 * @brief Some common ways to override #Callback::grp_name and #Callback::prt_name
 */

#ifndef CALLBACK_UTILS_NAME_HPP
#define CALLBACK_UTILS_NAME_HPP

#include <string>

#include "callback.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::grp_name and #Callback::prt_name
 */
namespace name
{

    /*! @brief Illustris-type hdf5 format.
     *
     * @tparam PartType     the particle type
     */
    template<typename AFields, uint8_t PartType>
    struct Illustris :
        virtual public Callback<AFields>
    {// {{{
        std::string grp_name () const override final
        {
            return "Group/";
        }
        std::string prt_name () const override final
        {
            return "PartType" + std::to_string(PartType)+"/";
        }
    };// }}}
    
    /*! @brief Illustris-type hdf5 format with custom rockstar.
     *
     * @tparam PartType     the particle type
     */
    template<typename AFields, uint8_t PartType>
    struct IllustrisRockstar :
        virtual public Callback<AFields>
    {// {{{
        std::string grp_name () const override final
        {
            return "";
        }
        std::string prt_name () const override final
        {
            return "PartType" + std::to_string(PartType)+"/";
        }
    };// }}}

    /*! @brief Gadget-type hdf5 format
     *
     * This is just #CallbackUtils::name::Illustris with PartType=1.
     */
    template<typename AFields>
    struct Gadget :
        virtual public Callback<AFields>,
        public Illustris<AFields, 1>
    { };

} // namespace name

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_NAME_HPP
