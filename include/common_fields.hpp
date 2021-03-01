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
} // namespace IllustrisFields

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
