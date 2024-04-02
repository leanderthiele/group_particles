/*! @file callback_utils_chunk.hpp
 *
 * @brief Some common ways to override #Callback::grp_chunk and #Callback::prt_chunk
 */

#ifndef CALLBACK_UTILS_CHUNK_HPP
#define CALLBACK_UTILS_CHUNK_HPP

#include <string>

#include "callback.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::grp_chunk and #Callback::prt_chunk
 */
namespace chunk {

    /*! @brief there is a single file containing all groups.
     */
    template<typename AFields>
    class SingleGrp :
        virtual public Callback<AFields>
    {// {{{
        const std::string grp_fname;
    public :
        /*! @param[in] grp_fname    the path to the single group file.
         */
        SingleGrp (const std::string &grp_fname_) :
            grp_fname(grp_fname_)
        { }

        bool grp_chunk (size_t chunk_idx, std::string &fname) const override final
        {
            fname = grp_fname;
            return !chunk_idx;
        }
    };// }}}

    /*! @brief there is a single file containing all particles.
     */
    template<typename AFields>
    class SinglePrt :
        virtual public Callback<AFields>
    {// {{{
        const std::string prt_fname;
    public :
        /*! @param[in] prt_fname    the path to the single particle file.
         */
        SinglePrt (const std::string &prt_fname_) :
            prt_fname(prt_fname_)
        { }
        bool prt_chunk (size_t chunk_idx, std::string &fname) const override final
        {
            fname = prt_fname;
            return !chunk_idx;
        }
    };// }}}

    /*! @brief there are multiple group chunks.
     */
    template<typename AFields>
    class MultiGrp :
        virtual public Callback<AFields>
    {// {{{
        const std::string grp_fname;
        const size_t min_idx, max_idx;
    public :
        /*! @param grp_fname    pattern to construct the group file names.
         *                      This should be a complete path with a placeholder (e.g. `"%lu"`)
         *                      that is replaced by the chunk_idx argument to #Callback::grp_chunk.
         *  @param max_idx      the last possible group chunk index (so there are max_idx+1 group chunks to read)
         *  @param min_idx      the frst possible group chunk index (zero by default)
         */
        MultiGrp (const std::string &grp_fname_, size_t max_idx_, size_t min_idx_=0UL) :
            grp_fname(grp_fname_), min_idx(min_idx_), max_idx(max_idx_)
        { }

        bool grp_chunk (size_t chunk_idx, std::string &fname) const override final
        {
            char buf[grp_fname.size()+10];
            std::sprintf(buf, grp_fname.c_str(), min_idx+chunk_idx);
            fname = std::string(buf);
            return min_idx+chunk_idx <= max_idx;
        }
    };// }}}

    /*! @brief there are multiple particle chunks.
     */
    template<typename AFields>
    class MultiPrt :
        virtual public Callback<AFields>
    {// {{{
        const std::string prt_fname;
        const size_t min_idx, max_idx;
    public :
        /*! @param prt_fname    pattern to construct the particle file names.
         *                      This should be a complete path with a placeholder (e.g. `"%lu"`)
         *                      that is replaced by the chunk_idx argument to #Callback::prt_chunk.
         *  @param max_idx      the last possible particle chunk index (so there are max_idx+1 particle chunks to read)
         *  @param min_idx      the first possible particle chunk index (zero by default)
         */
        MultiPrt (const std::string &prt_fname_, size_t max_idx_, size_t min_idx_=0UL) :
            prt_fname(prt_fname_), min_idx(min_idx_), max_idx(max_idx_)
        { }

        bool prt_chunk (size_t chunk_idx, std::string &fname) const override final
        {
            char buf[prt_fname.size()+10];
            std::sprintf(buf, prt_fname.c_str(), min_idx+chunk_idx);
            fname = std::string(buf);
            return min_idx+chunk_idx <= max_idx;
        }
    };// }}}

    /*! @brief wrapper for the common case that there is one group chunk and one particle chunk
     */
    template<typename AFields>
    struct Single :
        virtual public Callback<AFields>,
        public SingleGrp<AFields>, public SinglePrt<AFields>
    {// {{{
        /*! see the constructors for #SingleGrp and #SinglePrt
         */
        Single (const std::string &grp_fname,
                const std::string &prt_fname) :
            SingleGrp<AFields>(grp_fname), SinglePrt<AFields>(prt_fname) { }
    };// }}}

    /*! @brief wrapper for the common case that there are multiple group chunks and multiple particle chunks
     */
    template<typename AFields>
    struct Multi :
        virtual public Callback<AFields>,
        public MultiGrp<AFields>, public MultiPrt<AFields>
    {// {{{
        /*! see the constructors for #MultiGrp and #MultiPrt
         */
        Multi (const std::string &grp_fname, size_t grp_max_idx,
               const std::string &prt_fname, size_t prt_max_idx) :
            MultiGrp<AFields>(grp_fname, grp_max_idx),
            MultiPrt<AFields>(prt_fname, prt_max_idx) { }
    };// }}}

} // namespace chunk

} // namespace CallbackUtils

#endif // CALLBACK_UTILS_CHUNK_HPP
