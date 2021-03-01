/*! @file callback_utils_grp_action.hpp
 *
 * @brief Some common ways to override #Callback::grp_action.
 */


#ifndef CALLBACK_UTILS_GRP_ACTION_HPP
#define CALLBACK_UTILS_GRP_ACTION_HPP

#include <cassert>
#include <functional>
#include <vector>
#include <type_traits>

#include "callback.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::grp_action.
 */
namespace grp_action {

    /*! @brief base class to have multiple actions performed in #Callback::grp_action
     *
     * The user should not inherit from this class directly
     * but rather through the #CallbackUtils::grp_action::MultiGrpAction class.
     */
    template<typename AFields>
    class MultiGrpActionBase :
        virtual public Callback<AFields>
    {// {{{
        using typename Callback<AFields>::GrpProperties;
        static constexpr size_t buf_size = 64UL;
        size_t N_actions = 0UL;
        std::pair<void *, std::function<void(void *, const GrpProperties &)>> grp_actions[buf_size];
    protected :
        void register_fct (void *obj, std::function<void(void *, const GrpProperties &)> fct)
        {
            grp_actions[N_actions++] = std::make_pair(obj, fct);
            assert(N_actions < buf_size);
        }
    public :
        void grp_action (const GrpProperties &grp) override final
        {
            for (size_t ii=0; ii != N_actions; ++ii)
                grp_actions[ii].second(grp_actions[ii].first, grp);
        }
    };// }}}

    /*! @brief interface class to add an action to be performed during #Callback::grp_action
     *
     * @tparam Child        the inheriting type
     *
     * @note Child has to make this class friend.
     *
     * See #CallbackUtils::grp_action::StoreGrpHomogeneous for an example.
     */
    template<typename AFields, typename Child>
    class MultiGrpAction :
        virtual public Callback<AFields>,
        virtual private MultiGrpActionBase<AFields>
    {// {{{
        using typename Callback<AFields>::GrpProperties;
        static void this_grp_action_static (void *obj, const GrpProperties &grp)
        {
            Child *p = (Child *)obj;
            p->this_grp_action(grp);
        }
    protected :
        /*! The user should override this function with the desired action that is to be performed
         *  as part of the call to #Callback::grp_action.
         */
        virtual void this_grp_action (const GrpProperties &grp) = 0;
    public :
        MultiGrpAction ()
        {
            MultiGrpActionBase<AFields>::register_fct(this, this_grp_action_static);
        }
    };// }}}

    /*! @brief implements the common case of storing a data item of identical type for each group.
     *
     * @tparam Tdata        the data type that is to be stored for each group.
     */
    template<typename AFields, typename Tdata>
    class StoreGrpHomogeneous :
        virtual public Callback<AFields>,
        private MultiGrpAction<AFields, StoreGrpHomogeneous<AFields, Tdata>>
    {// {{{
        friend MultiGrpAction<AFields, StoreGrpHomogeneous<AFields, Tdata>>;
        using typename Callback<AFields>::GrpProperties;
        std::vector<Tdata> &data;
        void this_grp_action (const GrpProperties &grp) override final
        {
            data.push_back(grp_reduce(grp));
        }
    protected :
        /*! The user should override this function to implement how the properties of a group
         *  should be reduced to the desired Tdata type.
         *
         *  For each group, the result of this function will then be appended to the data vector.
         */
        virtual Tdata grp_reduce (const GrpProperties &grp) const = 0;
    public :
        /*! @param data     a zero-length vector that will be filled during the #Callback::grp_action calls
         *                  (according to the #grp_reduce result)
         */
        StoreGrpHomogeneous (std::vector<Tdata> &data_) :
            data(data_)
        {
            assert(data.empty());
        }
    };// }}}

    /*! @brief implements the common case of storing a single group property for each group.
     *
     * @tparam Field        the group property that is stored for each group
     * @tparam storeasT     the type the property should be cast to
     *
     * @note Currently, only 1-dimensional Fields are supported (and storeasT must be arithmetic)
     */
    template<typename AFields, typename Field, typename storeasT = typename Field::value_type>
    class StoreGrpProperty :
        virtual public Callback<AFields>,
        private StoreGrpHomogeneous<AFields, storeasT>
    {// {{{
        static_assert(Field::dim == 1, "Currently not implemented, could probably be done");
        static_assert(std::is_arithmetic_v<storeasT>);
        storeasT grp_reduce (const typename Callback<AFields>::GrpProperties &grp) const override final
        {
            return grp.template get<Field>();
        }
    public :
        /*! @param data     a zero-length vector that will be filled with the desired group property Field.
         */
        StoreGrpProperty (std::vector<storeasT> &data) :
            StoreGrpHomogeneous<AFields, storeasT>(data)
        { }
    };// }}}

} // namespace grp_action

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_GRP_ACTION_HPP
