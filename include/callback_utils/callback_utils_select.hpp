/*! @file callback_utils_select.hpp
 *
 * @brief Some common ways to override #Callback::grp_select
 */

#ifndef CALLBACK_UTILS_SELECT_HPP
#define CALLBACK_UTILS_SELECT_HPP

#include <cassert>
#include <utility>
#include <functional>

#include "callback.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::grp_select
 */
namespace select
{

    /*! @brief base class to apply multiple selection functions in #Callback::grp_select.
     *
     * The user should not inherit directly from this class
     * but rather through the #CallbackUtils::select::MultiSelect class.
     */
    template<typename AFields>
    class MultiSelectBase : virtual public Callback<AFields>
    {// {{{
        typedef typename Callback<AFields>::GrpProperties GrpProperties;
        static constexpr size_t buf_size = 64UL;
        size_t N_selects = 0UL;
        std::pair<void *, std::function<bool(void *, const GrpProperties &)>> selectors[buf_size];
    protected :
        void register_select (void *obj, std::function<bool(void *, const GrpProperties &)> fct)
        {
            selectors[N_selects++] = std::make_pair(obj, fct);
            assert(N_selects < buf_size);
        }
    public :
        bool grp_select (const GrpProperties &grp) const override final
        {
            for (size_t ii=0; ii != N_selects; ++ii)
                if (!selectors[ii].second(selectors[ii].first, grp))
                    return false;
            return true;
        }
    };// }}}

    /*! @brief interface class to include a selection function in #Callback::grp_select.
     *
     * @tparam Child        the inheriting type
     *
     * See #CallbackUtils::select::Window for an example.
     */
    template<typename AFields, typename Child>
    class MultiSelect : virtual public Callback<AFields>,
                        virtual private MultiSelectBase<AFields>
    {// {{{
        typedef typename Callback<AFields>::GrpProperties GrpProperties;
        static bool this_grp_select_static (void *obj, const GrpProperties &grp)
        {
            Child *p = (Child *)obj;
            return p->this_grp_select(grp);
        }
    protected :
        /*! The user should override this function with the desired selection that is to be
         *  performed on the group's properties.
         */
        virtual bool this_grp_select (const GrpProperties &grp) const = 0;
    public :
        MultiSelect ()
        {
            MultiSelectBase<AFields>::register_select(this, this_grp_select_static);
        }
    };// }}}

    /*! @brief select only groups that have some 1-dimensional property fall into an interval.
     *
     * @tparam Field    the group property that should be checked.
     */
    template<typename AFields, typename Field>
    class Window : virtual public Callback<AFields>,
                   public MultiSelect<AFields, Window<AFields, Field>>
    {// {{{
        typedef typename Callback<AFields>::GrpProperties GrpProperties;
        static_assert(Field::dim == 1);
        static_assert(Field::type == FieldTypes::GrpFld);
        static_assert(std::is_floating_point_v<typename Field::value_type>);
        typename Field::value_type min_val, max_val;
        bool this_grp_select (const GrpProperties &grp) const override final
        {
            auto x = grp.template get<Field>();
            return x > min_val && x < max_val;
        }
    public :
        /*! @param min_val      lower edge of the interval
         *  @param max_val      upper edge of the interval
         */
        Window (typename Field::value_type min_val_, typename Field::value_type max_val_) :
            min_val { min_val_ }, max_val  { max_val_ }
        { }
    };// }}}
    
    /*! @brief select only groups that have some 1-dimensional property above a certain value.
     *
     * @tparam Field    the group property that should be checked.
     */
    template<typename AFields, typename Field>
    struct LowCutoff : virtual public Callback<AFields>,
                       public Window<AFields, Field>
    {// {{{
        /*! @param min_val      lower limit on Field.
         */
        LowCutoff (typename Field::value_type min_val) :
            Window<AFields, Field> { min_val, std::numeric_limits<typename Field::value_type>::max() } { }
    };// }}}
    
    /*! @brief select only groups that have some 1-dimensional property below a certain value.
     *
     * @tparam Field    the group property that should be checked.
     */
    template<typename AFields, typename Field>
    struct HighCutoff : virtual public Callback<AFields>,
                        public Window<AFields, Field>
    {// {{{
        /*! @param max_val      upper limit on Field.
         */
        HighCutoff (typename Field::value_type max_val) :
            Window<AFields, Field> { std::numeric_limits<typename Field::value_type>::min(), max_val } { }
    };// }}}

} // namespace select

} // namespace CallbackUtils


#endif // CALLBACK_UTILS_SELECT_HPP
