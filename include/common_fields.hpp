/*! @file common_fields.hpp
 *
 * @brief Contains some common field types that can be used as template
 *        arguments for the #GrpFields and #PrtFields classes and for the
 *        #Callback::GrpProperties::get and #Callback::PrtProperties::get
 *        methods.
 *
 * They are constructed using the #FIELD macro defined in fields.hpp.
 * All types in this file are non-constructible and can only be used
 * as template arguments.
 */

#ifndef COMMON_FIELDS_HPP
#define COMMON_FIELDS_HPP

#include "fields.hpp"

/*! @brief Some fields contained in the Illustris simulation data products.
 */
namespace IllustrisFields
{
    // Particle fields
    FIELD(CenterOfMass, 3, float, FieldTypes::PrtFld, true);
    FIELD(Coordinates, 3, double, FieldTypes::PrtFld, true);
    FIELD(Velocities, 3, float, FieldTypes::PrtFld, false);
    FIELD(Density, 1, float, FieldTypes::PrtFld, false);
    FIELD(ElectronAbundance, 1, float, FieldTypes::PrtFld, false);
    FIELD(Masses, 1, float, FieldTypes::PrtFld, false);
    FIELD(InternalEnergy, 1, float, FieldTypes::PrtFld, false);
    FIELD(Potential, 1, float, FieldTypes::PrtFld, false);
    FIELD(StarFormationRate, 1, float, FieldTypes::PrtFld, false);
    FIELD(BH_CumEgyInjection_QM, 1, float, FieldTypes::PrtFld, false);
    FIELD(BH_CumEgyInjection_RM, 1, float, FieldTypes::PrtFld, false);
    // TODO

    // Group fields
    FIELD(GroupCM, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupPos, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupMass, 1, float, FieldTypes::GrpFld, false);
    FIELD(GroupVel, 3, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Crit200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Crit200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Crit500, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Crit500, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Mean200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Mean200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_TopHat200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_TopHat200, 1, float, FieldTypes::GrpFld, false);
    // TODO
} // namespace IllustrisFields

/*! @brief Fields contained in my hacky Rockstar data products -- not for public use.
 */
namespace RockstarFields
{
    FIELD(pos, 3, float, FieldTypes::GrpFld, true);
    FIELD(M200c, 1, float, FieldTypes::GrpFld, false);
    FIELD(R200c, 1, float, FieldTypes::GrpFld, false);
    FIELD(ang_mom, 3, float, FieldTypes::GrpFld, false);
    FIELD(Xoff, 1, float, FieldTypes::GrpFld, false);
    FIELD(Voff, 1, float, FieldTypes::GrpFld, false);
    FIELD(Vmax, 1, float, FieldTypes::GrpFld, false);
    FIELD(Vrms, 1, float, FieldTypes::GrpFld, false);
    FIELD(Rs, 1, float, FieldTypes::GrpFld, false);
    FIELD(vel, 3, float, FieldTypes::GrpFld, false);
    FIELD(Spin, 1, float, FieldTypes::GrpFld, false);
    FIELD(rs_klypin, 1, float, FieldTypes::GrpFld, false);
    FIELD(M200c_all, 1, float, FieldTypes::GrpFld, false);
    FIELD(Mvir, 1, float, FieldTypes::GrpFld, false);
    FIELD(M200b, 1, float, FieldTypes::GrpFld, false);
    FIELD(M500c, 1, float, FieldTypes::GrpFld, false);
    FIELD(M2500c, 1, float, FieldTypes::GrpFld, false);
    FIELD(spin_bullock, 1, float, FieldTypes::GrpFld, false);
    FIELD(b_to_a, 1, float, FieldTypes::GrpFld, false);
    FIELD(c_to_a, 1, float, FieldTypes::GrpFld, false);
} // namespace RockstarFields

/*! @brief Fields contained in my hacky Rockstar data products -- not for public use.
 */
namespace CustomRockstarFields
{
    FIELD(Pos, 3, float, FieldTypes::GrpFld, true);
    FIELD(Mvir, 1, float, FieldTypes::GrpFld, false);
    FIELD(ID, 1, int64_t, FieldTypes::GrpFld, false);
    FIELD(M_all, 1, float, FieldTypes::GrpFld, false);
    FIELD(ALT_M_0, 1, float, FieldTypes::GrpFld, false);
    FIELD(ALT_M_1, 1, float, FieldTypes::GrpFld, false);
    FIELD(ALT_M_2, 1, float, FieldTypes::GrpFld, false);
    FIELD(ALT_M_3, 1, float, FieldTypes::GrpFld, false);
    FIELD(Rvir, 1, float, FieldTypes::GrpFld, false);
    FIELD(J, 3, float, FieldTypes::GrpFld, false);
    FIELD(A, 3, float, FieldTypes::GrpFld, false);
    FIELD(A2, 3, float, FieldTypes::GrpFld, false);
    FIELD(X_off, 1, float, FieldTypes::GrpFld, false);
    FIELD(V_off, 1, float, FieldTypes::GrpFld, false);
    FIELD(Vmax, 1, float, FieldTypes::GrpFld, false);
    FIELD(Vrms, 1, float, FieldTypes::GrpFld, false);
    FIELD(Rs, 1, float, FieldTypes::GrpFld, false);
    FIELD(Vel, 3, float, FieldTypes::GrpFld, false);
    FIELD(Spin, 1, float, FieldTypes::GrpFld, false);
    FIELD(Bullock_Spin, 1, float, FieldTypes::GrpFld, false);
    FIELD(Klypin_Rs, 1, float, FieldTypes::GrpFld, false);
    FIELD(b_to_a, 1, float, FieldTypes::GrpFld, false);
    FIELD(c_to_a, 1, float, FieldTypes::GrpFld, false);
    FIELD(b_to_a2, 1, float, FieldTypes::GrpFld, false);
    FIELD(c_to_a2, 1, float, FieldTypes::GrpFld, false);
} // namespace CustomRockstarFields
  //
/*! @brief Fields contained in my Uchuu rockstar catalog
 */
namespace UchuuRockstarFields
{
    FIELD(Pos, 3, double, FieldTypes::GrpFld, true);
    FIELD(Mvir, 1, double, FieldTypes::GrpFld, false);
    FIELD(id, 1, int64_t, FieldTypes::GrpFld, false);
    FIELD(pid, 1, int64_t, FieldTypes::GrpFld, false);
    FIELD(M_vir_all, 1, double, FieldTypes::GrpFld, false);
    FIELD(M200b, 1, double, FieldTypes::GrpFld, false);
    FIELD(M200c, 1, double, FieldTypes::GrpFld, false);
    FIELD(M2500c, 1, double, FieldTypes::GrpFld, false);
    FIELD(M500c, 1, double, FieldTypes::GrpFld, false);
    FIELD(Rvir, 1, double, FieldTypes::GrpFld, false);
    // FIELD(J, 3, float, FieldTypes::GrpFld, false);
    // FIELD(A, 3, float, FieldTypes::GrpFld, false);
    // FIELD(A2, 3, float, FieldTypes::GrpFld, false);
    FIELD(X_off, 1, double, FieldTypes::GrpFld, false);
    FIELD(V_off, 1, double, FieldTypes::GrpFld, false);
    FIELD(vmax, 1, double, FieldTypes::GrpFld, false);
    FIELD(vrms, 1, double, FieldTypes::GrpFld, false);
    FIELD(rs, 1, double, FieldTypes::GrpFld, false);
    FIELD(Vel, 3, double, FieldTypes::GrpFld, false);
    FIELD(Spin, 1, double, FieldTypes::GrpFld, false);
    FIELD(Spin_Bullock, 1, double, FieldTypes::GrpFld, false);
    FIELD(Rs_Klypin, 1, double, FieldTypes::GrpFld, false);
    FIELD(b_to_a, 1, double, FieldTypes::GrpFld, false);
    FIELD(c_to_a, 1, double, FieldTypes::GrpFld, false);
    FIELD(b_to_a2, 1, double, FieldTypes::GrpFld, false);
    FIELD(c_to_a2, 1, double, FieldTypes::GrpFld, false);
} // namespace UchuuRockstarFields

/*! @brief Some fields contained in the SIMBA simulation data products.
 *
 * These are basically the same as the IllustrisFields with the difference
 * that the Coordinates are stored as 32 bit floats in SIMBA
 * (same as in the Illustris mini snapshots)
 */
namespace SIMBAFields
{
    // Particle fields
    FIELD(CenterOfMass, 3, float, FieldTypes::PrtFld, true);
    FIELD(Coordinates, 3, float, FieldTypes::PrtFld, true);
    FIELD(Density, 1, float, FieldTypes::PrtFld, false);
    FIELD(ElectronAbundance, 1, float, FieldTypes::PrtFld, false);
    FIELD(Masses, 1, float, FieldTypes::PrtFld, false);
    FIELD(InternalEnergy, 1, float, FieldTypes::PrtFld, false);
    // TODO

    // Group fields
    FIELD(GroupCM, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupPos, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupMass, 1, float, FieldTypes::GrpFld, false);
    FIELD(GroupVel, 3, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Crit200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Crit200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Crit500, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Crit500, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_Mean200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_Mean200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_R_TopHat200, 1, float, FieldTypes::GrpFld, false);
    FIELD(Group_M_TopHat200, 1, float, FieldTypes::GrpFld, false);
    // TODO
} // namespace SIMBAFields

/*! @brief Fields contained in a typical Gadget simulation.
 *
 * Care should be taken with the value types, depending on the Gadget
 * build these can either be double or float.
 */
namespace GadgetFields
{
    // Particle fields
    FIELD(Coordinates, 3, float, FieldTypes::PrtFld, true);
    FIELD(Velocities, 3, float, FieldTypes::PrtFld, false);

    // Group fields
    FIELD(GroupPos, 3, float, FieldTypes::GrpFld, true);
    FIELD(GroupVel, 3, float, FieldTypes::GrpFld, false);
    FIELD(GroupMass, 1, float, FieldTypes::GrpFld, false);
    FIELD(GroupAscale, 1, float, FieldTypes::GrpFld, false);
} // namespace GadgetFields



#endif // COMMON_FIELDS_HPP
