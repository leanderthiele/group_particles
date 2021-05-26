/*! @file callback_utils_meta_init.hpp
 *
 * @brief Some common ways to override #Callback::read_grp_meta_init and #Callback::read_prt_meta_init
 */


#ifndef CALLBACK_UTILS_META_INIT_HPP
#define CALLBACK_UTILS_META_INIT_HPP

#include <cassert>
#include <functional>
#include <memory>

#include "callback.hpp"
#include "hdf5_utils.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::read_grp_meta_init and #Callback::read_prt_meta_init
 */
namespace meta_init
{

    /*! @brief base class to have multiple actions performed in #Callback::read_prt_meta_init
     *
     * The user should not inherit directly from this class
     * but rather through the #CallbackUtils::meta_init::MultiPrtMetaInit class.
     */
    template<typename AFields>
    class MultiPrtMetaInitBase :
        virtual public Callback<AFields>
    {// {{{
        static constexpr size_t buf_size = 64UL;
        size_t N_inits = 0UL;
        std::pair<void *, std::function<void(void *, std::shared_ptr<H5::H5File>)>> inits[buf_size];
    protected :
        void register_init (void *obj, std::function<void(void *, std::shared_ptr<H5::H5File>)> fct)
        {
            inits[N_inits++] = std::make_pair(obj, fct);
            assert(N_inits < buf_size);
        }
    public :
        void read_prt_meta_init (std::shared_ptr<H5::H5File> fptr) override final
        {
            for (size_t ii=0; ii != N_inits; ++ii)
                inits[ii].second(inits[ii].first, fptr);
        }
    };// }}}

    /*! @brief interface class to add an action to be performed during #Callback::read_prt_meta_init
     *
     * @tparam Child        the inheriting type
     *
     * @note Child has to make this class friend.
     *
     * See #CallbackUtils::meta_init::IllustrisCosmology for an example.
     */
    template<typename AFields, typename Child>
    class MultiPrtMetaInit :
        virtual public Callback<AFields>,
        virtual private MultiPrtMetaInitBase<AFields>
    {// {{{
        static void this_prt_meta_init_static (void *obj, std::shared_ptr<H5::H5File> fptr)
        {
            Child *p = (Child *)obj;
            p->this_prt_meta_init(fptr);
        }
    protected :
        /*! The user should override this function with the desired action that is to be performed
         *  as part of the call to #Callback::read_prt_meta_init.
         */
        virtual void this_prt_meta_init (std::shared_ptr<H5::H5File> fptr) = 0;
    public :
        MultiPrtMetaInit ()
        {
            MultiPrtMetaInitBase<AFields>::register_init(this, this_prt_meta_init_static);
        }
    };// }}}

    /*! @brief stores some cosmology-related meta-data from an Illustris-type header.
     */
    template<typename AFields>
    class IllustrisCosmology :
        virtual public Callback<AFields>,
        private MultiPrtMetaInit<AFields, IllustrisCosmology<AFields>>
    {// {{{
        friend MultiPrtMetaInit<AFields, IllustrisCosmology<AFields>>;
        void this_prt_meta_init (std::shared_ptr<H5::H5File> fptr) override final
        {
            auto header = fptr->openGroup("/Header");
            #define READ(x) x = hdf5Utils::read_scalar_attr<double,double>(header, #x)
            READ(HubbleParam);
            READ(Omega0);
            READ(OmegaLambda);
            READ(OmegaBaryon);
            READ(Redshift);
            READ(Time);
            #undef READ
        }
    protected :
        double HubbleParam, /*!< @brief Hubble parameter */
               Omega0, /*!< @brief matter density */
               OmegaLambda, /*!< @brief dark energy density */
               OmegaBaryon, /*!< @brief baryonic density */
               Redshift, /*!< @brief redshift */
               Time; /*!< @brief scale factor */
    };// }}}

    /*! @brief stores the MassTable from an Illustris-type header.
     */
    template<typename AFields>
    class IllustrisMassTable :
        virtual public Callback<AFields>,
        private MultiPrtMetaInit<AFields, IllustrisMassTable<AFields>>
    {// {{{
        friend MultiPrtMetaInit<AFields, IllustrisCosmology<AFields>>;
        void this_prt_meta_init (std::shared_ptr<H5::H5File> fptr) override final
        {
            auto header = fptr->openGroup("/Header");
            Ntypes = hdf5Utils::read_vector_attr<double,double>(header, "MassTable", MassTable);
        }
    protected :
        double MassTable[16]; /*!< @brief the mass table (length #Ntypes), indices corresponding to particle types */
        size_t Ntypes; /*!< @brief length of the mass table */
    };// }}}

} // namespace meta_init


} // namespace CallbackUtils


#endif // CALLBACK_UTILS_META_INIT_HPP
