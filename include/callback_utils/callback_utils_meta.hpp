/*! @file callback_utils_meta.hpp
 *
 * @brief Some common ways to override #Callback::read_grp_meta and #Callback::read_prt_meta
 */


#ifndef CALLBACK_UTILS_META_HPP
#define CALLBACK_UTILS_META_HPP

#include "callback.hpp"
#include "hdf5_utils.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::read_grp_meta and #Callback::read_prt_meta
 */
namespace meta
{

    /*! @brief retrieves meta-data from an Illustris-type simulation.
     *
     * @tparam PartType     the particle type
     */
    template<typename AFields, uint8_t PartType>
    struct Illustris :
        virtual public Callback<AFields>
    {// {{{
        void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                            size_t &Ngroups) const override final
        {
            auto header = fptr->openGroup("/Header");
            Ngroups = hdf5Utils::read_scalar_attr<int32_t,size_t>(header, "Ngroups_ThisFile");
            header.close();
        }

        void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                            coord_t &Bsize, size_t &Npart) const override final
        {
            auto header = fptr->openGroup("/Header");
            Bsize = hdf5Utils::read_scalar_attr<double,coord_t>(header, "BoxSize");
            Npart = hdf5Utils::read_vector_attr<int32_t,size_t>(header, "NumPart_ThisFile", PartType);
            header.close();
        }
    };// }}}

    /*! @brief retrieves meta-data from an Illustris-type simulation with my custom rockstar hdf5.
     *
     * @tparam PartType     the particle type
     */
    template<typename AFields, uint8_t PartType>
    struct IllustrisRockstar :
        virtual public Callback<AFields>
    {// {{{
        void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                            size_t &Ngroups) const override final
        {
            auto header = fptr->openGroup("/Header");
            Ngroups = hdf5Utils::read_scalar_attr<int64_t,size_t>(header, "Nhosts");
            header.close();
        }

        void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                            coord_t &Bsize, size_t &Npart) const override final
        {
            auto header = fptr->openGroup("/Header");
            Bsize = hdf5Utils::read_scalar_attr<double,coord_t>(header, "BoxSize");
            Npart = hdf5Utils::read_vector_attr<int32_t,size_t>(header, "NumPart_ThisFile", PartType);
            header.close();
        }
    };// }}}

    /*! @brief retrieves meta-data from a SIMBA-type simulation.
     *
     * This is essentially the same as for an Illustris type sim (i.e. #CallbackUtils::meta::Illustris),
     * but with the subtle difference that NumPart_ThisFile is stored as a uint32_t
     *
     * @tparam PartType     the particle type
     */
    template<typename AFields, uint8_t PartType>
    struct SIMBA :
        virtual public Callback<AFields>
    {// {{{
        void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                            size_t &Ngroups) const override final
        {
            auto header = fptr->openGroup("/Header");
            Ngroups = hdf5Utils::read_scalar_attr<int32_t,size_t>(header, "Ngroups_ThisFile");
            header.close();
        }

        void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                            coord_t &Bsize, size_t &Npart) const override final
        {
            auto header = fptr->openGroup("/Header");
            Bsize = hdf5Utils::read_scalar_attr<double,coord_t>(header, "BoxSize");
            Npart = hdf5Utils::read_vector_attr<uint32_t,size_t>(header, "NumPart_ThisFile", PartType);
            header.close();
        }
    };// }}}

    /*! @brief retrieves meta-data from a Gadget-type simulation.
     *
     * This is just #CallbackUtils::meta::Illustris with PartType=1.
     */
    template<typename AFields>
    struct Gadget :
        virtual public Callback<AFields>,
        public Illustris<AFields, 1>
    { };

} // namespace meta

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_META_HPP
