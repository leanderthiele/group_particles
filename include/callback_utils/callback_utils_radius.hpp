/*! @file callback_utils_radius.hpp
 *
 * @brief Some common ways to override #Callback::grp_radius
 */


#ifndef CALLBACK_UTILS_RADIUS_HPP
#define CALLBACK_UTILS_RADIUS_HPP

#include "callback.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::grp_radius
 */
namespace radius
{

    /*! @brief radius is propertional to one of the group properties.
     *
     * @tparam RField       the group property the radius is proportional to
     */
    template<typename AFields, typename RField>
    class Simple : virtual public Callback<AFields>
    {// {{{
        static_assert(RField::dim == 1);
        static_assert(RField::type == FieldTypes::GrpFld);
        static_assert(std::is_floating_point_v<typename RField::value_type>);
        typename RField::value_type scaling;
    public :
        /*! @param scaling      the propertionality factor, the group radius is computed as scaling * RField
         */
        Simple (typename RField::value_type scaling_) :
            scaling { scaling_ }
        { }

        /*! default constructor sets the proportionality factor to 1.
         */
        Simple () :
            scaling { (typename RField::value_type)1.0 }
        { }

        coord_t grp_radius (const typename Callback<AFields>::GrpProperties &grp) const override final
        {
            return scaling * grp.template get<RField>();
        }
    };// }}}

} // namespace radius

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_RADIUS_HPP
