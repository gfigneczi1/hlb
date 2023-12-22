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

#ifndef DC_EMG_PDPLANECHANGE_HPP_INCLUDED
#define DC_EMG_PDPLANECHANGE_HPP_INCLUDED

namespace Rb
{
namespace Vmc
{

/// calculate level of dynamics of the given driver if sets parameters accordingly

class pdpLaneChange
{

public:
   /// This is the main method to optimize the target path
   ///
   /// The method will call all required methods to configure and solve the optimization problem.
   /// During the optimization the vector of paths m_pathsToEvaluate will iteratively (in a loop) be modified.
   /// In each iteration of the loop all vector elements will be evaluated and calculated, in case the are not
   /// available from a previous iteration. This means there will be roughly (number of iterations)*(vector
   /// elements) paths calculated and evaluated.
   ///
   /// \param [in] vehicle state including yawrate, velocity, road wheel angle
   /// \param [in] Frenet-frame path coefficients
   /// \param [in] parameter class of TrpLar

   void run(
      float input);

private:

   float                  m_variable{0.0F};
};  // end of class pdpLaneChange


}  // end of namespace Emg
}  // end of namespace Dc


#endif  // DC_EMG_MULTIPLESECTIONOPTIMIZER_HPP_INCLUDED
