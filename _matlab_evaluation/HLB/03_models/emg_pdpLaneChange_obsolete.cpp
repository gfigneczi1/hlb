///=============================================================================
///  C O P Y R I G H T
///-----------------------------------------------------------------------------
/// @copyright (c) 2022 by Robert Bosch GmbH. All rights reserved.
///
///  The reproduction, distribution and utilization of this file as
///  well as the communication of its contents to others without express
///  authorization is prohibited. Offenders will be held liable for the
///  payment of damages. All rights reserved in the event of the grant
///  of a patent, utility model or design.
///  @file
///=============================================================================

#include "emg_pdpLaneChange.hpp"

namespace Rb
{
namespace Vmc
{


void pdpLaneChange::run(
   float input)
{
   m_variable = input*2.0f;
}

}  // namespace Vmc
}  // namespace Rb
