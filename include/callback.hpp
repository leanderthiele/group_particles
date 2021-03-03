/*! @file callback.hpp
 *
 * @brief Contains the abstract base class that the user should subclass from
 *        in order to define the desired functionality of the #group_particles function.
 *
 */

#ifndef CALLBACK_HPP
#define CALLBACK_HPP

#include <cstddef>
#include <string>
#include <memory>
#include <type_traits>

#include "H5Cpp.h"

#include "fields.hpp"

/*! @brief The abstract base class the user should inherit from.
 * 
 * @tparam AFields a type constructed from the #AllFields template.
 *
 * By overriding this classes methods, the user defines the functionality of the code.
 *
 * See the documentation for #group_particles for the order in which the non-const
 * methods are called.
 */
template<typename AFields>
struct Callback
{
    /*! @brief Type describing either a group or a particle.
     *
     * @tparam T    either AFields::GroupFields or AllFields::ParticleFields
     *
     * This type contains all the information about a group or particle
     * read from the data files (with the AFields template parameter controlling
     * which fields are to be read).
     * User should only use the #get method to retrieve the information
     * they need.
     *
     * Specializations are in #GrpProperties and #PrtProperties.
     */
    template<typename T>
    class BaseProperties
    {// {{{
        // declare as char so we can advance byte-wise
        // (char is guaranteed to be one byte wide)
        char *data[T::Nfields];

        template<typename Field>
        static constexpr size_t get_field_idx ();

    public :
        // optional offset initializes at offset particles inside the data_in_memory pointers
        BaseProperties (void *data_in_memory[T::Nfields], size_t offset=0UL);

        // advances the internal pointers by the respective strides
        // User should not call this!
        void advance ();

        // retrieves the internal pointers
        // User should not need to call this
        void *operator[] (size_t idx) const;

        /*! Special case of the #get method, returning the coordinate field.
         */
        auto coord () const;
        
        /*! The only member function the user should need to call.
         *
         * @tparam Field    the field whose value should be retrieved.
         *
         * @return If the Field is 1-dimensional (e.g. a group mass), the value will be returned.
         *         Otherwise (e.g. for a particle velocity), a pointer to the first element will
         *         be returned.
         */
        template<typename Field>
        auto get () const;
    };// }}}

    /*! @brief Specialization of the #Callback::BaseProperties type to groups. */
    using GrpProperties = BaseProperties<typename AFields::GroupFields>;

    /*! @brief Specialization of the #Callback::BaseProperties type to particles. */
    using PrtProperties = BaseProperties<typename AFields::ParticleFields>;

    /*! @brief Where to find the group files.
     *
     * @param[in] chunk_idx     index of the group chunk the code wants to read (starting from 0).
     * @param[out] fname        file name corresponding to this index.
     *
     * @return true if there is a chunk corresponding to chunk_idx, false otherwise.
     *         By definition, must return true for chunk_idx=0.
     *
     * @note see #CallbackUtils::chunk for some overrides.
     */
    virtual bool grp_chunk (size_t chunk_idx, std::string &fname) const = 0;

    /*! @brief Where to find the particle files.
     *
     * @param[in] chunk_idx     index of the particle chunk the code wants to read (starting from 0).
     * @param[out] fname        file name corresponding to this index.
     *
     * @return true if there is a chunk corresponding to chunk_idx, false otherwise.
     *         By definition, must return true for chunk_idx=0.
     *
     * @note see #CallbackUtils::chunk for some overrides.
     */
    virtual bool prt_chunk (size_t chunk_idx, std::string &fname) const = 0;

    /*! @brief Where to find the group fields in the hdf5 file.
     *
     * @remark for Illustris-type formats, "Group/"
     *
     * @note see #CallbackUtils::name for some overrides.
     */
    virtual std::string grp_name () const = 0;

    /*! @brief Where to find the particle fields in the hdf5 file.
     *
     * @remark for Illustris-type formats, "PartType.../"
     *
     * @note see #CallbackUtils::name for some overrides.
     */
    virtual std::string prt_name () const = 0;

    /*! @brief Allows the user to read meta-data from the 0th group chunk.
     *
     * @param[in] fptr      Points to the opened 0th group chunk.
     *
     * @remark This function is trivially implemented, so does not need to be overriden.
     *
     * @note see #CallbackUtils::meta_init for some overrides.
     */
    virtual void read_grp_meta_init (std::shared_ptr<H5::H5File> fptr) { return; }

    /*! @brief Allows the user to read meta-data from the 0th particle chunk.
     *
     * @param[in] fptr      Points to the opened 0th particle chunk.
     *
     * @remark This function is trivially implemented, so does not need to be overriden.
     *
     * @note see #CallbackUtils::meta_init for some overrides.
     */
    virtual void read_prt_meta_init (std::shared_ptr<H5::H5File> fptr) { return; }

    /*! @brief Inform the code how many groups there are in a group chunk.
     *
     * @param[in] chunk_idx     index of the group chunk, starting from 0.
     * @param[in] fptr          pointer to the opened group chunk file.
     * @param[out] Ngroups      to be filled with the number of groups in this file.
     *
     * @remark In principle, the code could infer Ngroups from the size of the data arrays
     *         stored in the file, but it is safer to have this additional check.
     *
     * @note see #CallbackUtils::meta for some overrides.
     */
    virtual void read_grp_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                                size_t &Ngroups) const = 0;

    /*! @brief Inform the code how large the box is and how many particles there are in a particle chunk.
     *
     * @param[in] chunk_idx     index of the particle chunk, starting from 0.
     * @param[in] fptr          pointer to the opened particle chunk file.
     * @param[out] Bsize        to be filled with the size of the simulation box.
     *                          This should always be the same;
     *                          the code checks this assumption as an additional safety feature.
     *                          Units are the same as the ones used in the data arrays.
     * @param[out] Nparts       to be filled with the number of particles in this file.
     *
     * @remark In principle, the code could infer Nparts from the size of the data arrays
     *         stored in the file, but it is safer to have this additional check.
     *
     * @note see #CallbackUtils::meta for some overrides.
     */
    virtual void read_prt_meta (size_t chunk_idx, std::shared_ptr<H5::H5File> fptr,
                                coord_t &Bsize, size_t &Nparts) const = 0;

    /*! @brief Inform the code whether a group should be considered.
     *
     * @param[in] grp       properties of this group.
     * 
     * @return whether this group should be considered.
     *
     * @remark This function is trivially implemented, so does not need to be overriden.
     *
     * @note see #CallbackUtils::select for some overrides.
     */
    virtual bool grp_select (const GrpProperties &grp) const { return true; }

    /*! @brief Action to take for each group for which #grp_select returned true.
     *
     * @param[in] grp       properties of this group.
     *
     * @remark This function will be called in the order of groups encountered,
     *         corresponding to the grp_idx argument in the #prt_action method.
     * @warning There can be a cross-reaction with certain overrides of the #prt_action
     *          method. In that case, overriding this method directly will lead to a
     *          compiler error. It is advised that the user always uses inheritance
     *          from the #CallbackUtils::grp_action::MultiGrpAction interface
     *          to indirectly override this method.
     *
     * @note see #CallbackUtils::grp_action for some overrides.
     */
    virtual void grp_action (const GrpProperties &grp) = 0;

    /*! @brief Inform the code how large this group is.
     *
     * @param[in] grp       properties of this group.
     *
     * @return the radius of this group. Only particles falling within
     *         this distance from the group coordinate will be passed
     *         to the #prt_action method.
     *
     * @note see #CallbackUtils::radius for some overrides.
     */
    virtual coord_t grp_radius (const GrpProperties &grp) const = 0;

    /*! @brief Action to take for each particle that falls within #grp_radius from
     *         a group.
     *
     *  @param[in] grp_idx      index of this group, corresponding to the order in which
     *                          #grp_action was called.
     *  @param[in] grp          properties of this group. Well-designed code should not
     *                          need to use this argument, as the required data products
     *                          can be reduced more efficiently in the #grp_action method.
     *  @param[in] prt          properties of this particle.
     *  @param[in] Rsq          squared distance between this particle's coordinate and this
     *                          group's coordinate. The code already computes this so the
     *                          user should be able to use it without re-computing it.
     *
     * @note see #CallbackUtils::prt_action for some overrides.
     */
    virtual void prt_action (size_t grp_idx, const GrpProperties &grp,
                             const PrtProperties &prt, coord_t Rsq) = 0;
};


// --- Implementation of the BaseProperties struct ---

template<typename AFields>
template<typename T>
Callback<AFields>::BaseProperties<T>::BaseProperties (void *data_in_memory[T::Nfields], size_t offset)
{
    for (size_t ii=0; ii != T::Nfields; ++ii)
        data[ii] = (char *)(data_in_memory[ii]) + offset * T::strides_fcoord[ii];
}

template<typename AFields>
template<typename T>
inline void
Callback<AFields>::BaseProperties<T>::advance ()
{
    for (size_t ii=0; ii != T::Nfields; ++ii)
        data[ii] += T::strides_fcoord[ii];
}

template<typename AFields>
template<typename T>
inline void *
Callback<AFields>::BaseProperties<T>::operator[] (size_t idx) const
{
    return data[idx];
}

template<typename AFields>
template<typename T>
inline auto
Callback<AFields>::BaseProperties<T>::coord () const
{
    return (coord_t *)data[0];
}

template<typename AFields>
template<typename T>
template<typename Field>
constexpr size_t
Callback<AFields>::BaseProperties<T>::get_field_idx ()
{
    static_assert(T::field_type == Field::type);

    if constexpr (T::field_type == FieldTypes::GrpFld)
        return AFields::GroupFields::template idx<Field>;
    else
        return AFields::ParticleFields::template idx<Field>;
}

template<typename AFields>
template<typename T>
template<typename Field>
inline auto
Callback<AFields>::BaseProperties<T>::get () const
{
    constexpr size_t idx = get_field_idx<Field> ();

    if constexpr (Field::dim == 1)
        return (typename Field::value_type) *(typename Field::value_type *)data[idx];
    else
        if constexpr (Field::coord)
        // need to use the global coord_t that this field has been converted to
            return (coord_t *)data[idx];
        else
            return (typename Field::value_type *)data[idx];
}

#endif // CALLBACK_HPP
