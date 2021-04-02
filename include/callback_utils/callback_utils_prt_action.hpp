/*! @file callback_utils_prt_action.hpp
 *
 * @brief Some common ways to override #Callback::prt_action
 */


#ifndef CALLBACK_UTILS_PRT_ACTION_HPP
#define CALLBACK_UTILS_PRT_ACTION_HPP

#include "callback.hpp"
#include "callback_utils_grp_action.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::prt_action
 */
namespace prt_action {

    /*! @brief implements the common case that each group gets an item of identical type
     *         into which the particles belonging to the group are inserted in some way
     *
     *  @tparam Tdata       the data type to be stored for each group.
     *                      If Tdata implements a constructor of the signature
     *                      `Tdata::Tdata (const GrpProperties &)`
     *                      this constructor will be called with the group each data item
     *                      corresponds to (during the calls to #Callback::grp_action).
     *                      This allows the user to store some group data that they require
     *                      to insert particles into the data item.
     *                      If this behaviour is desired, the default constructor of Tdata
     *                      should be deleted (this is simply an extra code check).
     *                      If no constructor from a `const GrpProperties &` is found,
     *                      the default constructor will be called.
     */
    template<typename AFields, typename Tdata>
    class StorePrtHomogeneous :
        virtual public Callback<AFields>,
        private grp_action::MultiGrpAction<AFields, StorePrtHomogeneous<AFields, Tdata>>
    {// {{{
        friend grp_action::MultiGrpAction<AFields, StorePrtHomogeneous<AFields, Tdata>>;
        using typename Callback<AFields>::GrpProperties;
        using typename Callback<AFields>::PrtProperties;

        // some functionality to figure out whether Tdata has the method
        // void prt_insert (size_t grp_idx, const GrpProperties &grp,
        //                  const PrtProperties &prt, coord_t Rsq)
        // from https://stackoverflow.com/questions/87372/check-if-a-class-has-a-member-function-of-a-given-signature
        template<typename, typename T>
        struct has_prt_insert { static_assert(std::integral_constant<T, false>::value); };

        template<typename C, typename Ret, typename... Args>
        class has_prt_insert<C, Ret(Args...)>
        {
            template<typename T>
            static constexpr auto check (T *)
                -> typename std::is_same<
                                decltype(std::declval<T>().prt_insert(std::declval<Args>()...)),
                                Ret
                            >::type;
            
            template<typename>
            static constexpr std::false_type check(...);

            using type = decltype(check<C>(0));
        public :
            static constexpr bool value = type::value;
        };


        std::vector<Tdata> &data;

        void this_grp_action (const GrpProperties &grp) override final
        {
            if constexpr (std::is_constructible_v<Tdata, const GrpProperties &>)
            {
                static_assert(!std::is_default_constructible_v<Tdata>);
                data.emplace_back(grp);
            }
            else
                data.emplace_back();
        }
    protected :
        /*! The user should override this function to implement how a particle should be inserted
         *  into a Tdata item in the data vector.
         *
         * @note Alternatively, in case `Tdata` is a composite type, you can implement a function
         *       @code
         *       void Tdata::prt_insert (size_t grp_idx, const GrpProperties &grp,
         *                               const PrtProperties &prt, coord_t Rsq);
         *       @endcode
         *       In that case, this method does not need to be overriden.
         *
         * @param[in] grp_idx           index of this group, corresponding to the order of
         *                              #Callback::grp_action calls.
         * @param[in] grp               properties of this group.
         * @param[in] prt               properties of this particle.
         * @param[in] Rsq               squared distance between the particle's coordinate
         *                              and the group's coordinate.
         * @param[in,out] data_item     element in the data vector corresponding to this group
         *                              that is to be modified.
         *
         * @note Well-designed code should not require the grp_idx and grp arguments.
         *       Rather, if some properties of the groups are required in order to insert the particles,
         *       `Tdata` should store them by implementing a custom constructor as mentioned in the class docs.
         *
         */
        virtual void prt_insert (size_t grp_idx, const GrpProperties &grp,
                                 const PrtProperties &prt, coord_t Rsq,
                                 Tdata &data_item)
        {
            assert(("Need to override CallbackUtils::StorePrtHomogeneous<AFields,Tdata>::prt_insert if Tdata does not implement prt_insert method.", false));
        }

    public :
        /*! @param data     a zero-length vector whose elements will be constructed
         *                  (depending on which constructors Tdata has, see class docs)
         *                  during the calls to #Callback::grp_action.
         *                  The elements will then be modified during the calls to
         *                  #Callback::prt_action according to the user's implementation
         *                  of the prt_insert function.
         */
        StorePrtHomogeneous (std::vector<Tdata> &data_) :
            data(data_)
        { }

        void prt_action (size_t grp_idx, const GrpProperties &grp,
                         const PrtProperties &prt, coord_t Rsq) override final
        {
            if constexpr (has_prt_insert<Tdata, void(size_t, const GrpProperties &,
                                                     const PrtProperties &, coord_t)>::value)
                data[grp_idx].prt_insert(grp_idx, grp, prt, Rsq);
            else
                prt_insert(grp_idx, grp, prt, Rsq, data[grp_idx]);
        }
    };// }}}

} // namespace prt_action


} // namespace CallbackUtils

#endif // CALLBACK_UTILS_PRT_ACTION_HPP
