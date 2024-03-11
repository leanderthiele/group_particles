/*! @file callback_utils_prt_modify.hpp
 *
 *  @brief Some common ways to override #Callback::prt_modify
 */

#ifndef CALLBACK_UTILS_PRT_MODIFY_HPP
#define CALLBACK_UTILS_PRT_MODIFY_HPP

#include "callback.hpp"

namespace CallbackUtils {

/*! @brief Some common ways to override #Callback::prt_modify
 */
namespace prt_modify {


    template<typename AFields, typename VField, bool sqrta>
    class PrtRSD :
        virtual public Callback<AFields>
    {
        using typename Callback<AFields>::PrtProperties;
        coord_t rsd_factor;
        size_t rsd_direction;
    
    public :
        PrtRSD () = delete;

        PrtRSD (char rsd_direction_, double Omega_m, double z) {
            rsd_direction = rsd_direction_ - 'x';
            rsd_factor = (1.0+z)
                         / ( 100.0 * std::sqrt( Omega_m*std::pow(1+z, 3) + (1.0-Omega_m) ) );
            if constexpr (sqrta) rsd_factor /= std::sqrt(1.0+z);
        }

        void prt_modify (PrtProperties &prt) override final {
            const auto v = prt.template get<VField>()[rsd_direction];
            auto x = prt.coord()[rsd_direction];
            x += rsd_factor * v;
            // make sure periodicity still respected
            if (x > prt.Bsize) x -= prt.Bsize;
            else if (x < 0) x += prt.Bsize;
            prt.coord()[rsd_direction] = x;
        }
    };
}

}



#endif // CALLBACK_UTILS_PRT_MODIFY_HPP

