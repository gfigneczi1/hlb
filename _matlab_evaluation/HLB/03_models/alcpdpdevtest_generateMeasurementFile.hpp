//============================================================================================================
// C O P Y R I G H T
//------------------------------------------------------------------------------------------------------------
/// \copyright (C) 2021 Robert Bosch GmbH.
//
// The reproduction, distribution and utilization of this file as
// well as the communication of its contents to others without express
// authorization is prohibited. Offenders will be held liable for the
// payment of damages. All rights reserved in the event of the grant
// of a patent, utility model or design.
//============================================================================================================
/*
 * Original functionality by Schmid Steffen (XC-DX/EYM5)
*/

#ifndef DC_ALCPDPDEVTEST_GENERATEMEASUREMENTFILE_HPP_INCLUDED
#define DC_ALCPDPDEVTEST_GENERATEMEASUREMENTFILE_HPP_INCLUDED

#include "dc_interfaces/dc_utils/normalization/dc_norms_convenienttypes.hpp"

#include <stdio.h>
#include <vector>

namespace Rb
{
namespace Vmc
{
namespace DevTest
{

//============================================================================================================
/// \brief This struct is used to define a measurement variable.
/// The measurement variable value has to be defined once.
/// Pointers are used to access the variables while the test/simulation is running
/// \tparam [in] T - typename (Float, Bool, Int) is required for the variable pointer
//------------------------------------------------------------------------------------------------------------
template <typename T>
struct MeasurementValue
{
   //============================================================================================================
   /// \brief Construct a new Measurement Value object
   /// \param [in] nameIn - Name of the measurement value in the measurement file
   /// \param [in] unitIn - Unit of the measurement value
   /// \param [in] valueAddressIn - Pointer to the measurement value
   /// \param [in] idIn - ID of the measurement value, which has to be in increasing order
   /// \param [in] enumspecIn - optional: Define enum specification
   //------------------------------------------------------------------------------------------------------------
   MeasurementValue(
      std::string nameIn, std::string unitIn, T* valueAddressIn, int idIn, std::string enumspecIn = "")
      : name(nameIn)
      , unit(unitIn)
      , valueAddress(valueAddressIn)
      , id(idIn)
      , enumspec(enumspecIn){};

   T*          valueAddress{nullptr};
   std::string name{""};
   std::string unit{""};
   int         id{0};
   std::string enumspec = "";
};

//============================================================================================================
/// \brief This class is used to generate a measurement file while a test/simulation is running.
/// The defineMeasurementValues methods shall be called once to create a list of variables that need to
/// be measured. It is also possible to directly push_back measurement values to the corresponding lists.
/// The writeDataD97 method has to be called each cycle to collect the measurement data
//------------------------------------------------------------------------------------------------------------
class GenerateMeasurementFile
{
public:
   //============================================================================================================
   /// \brief Method to write the d97 file
   /// \param [in] cycleTime - CycleTime of the test - defines the time vector
   //------------------------------------------------------------------------------------------------------------
   void writeDataD97(const vfc::CSI::si_second_f32_t& cycleTime);
	void GenerateMeasurementFile::closeFile();

   //============================================================================================================
   /// \brief Method to generate the measurement file
   /// \param [in] testName   - Testname which also defines the output file name
   //------------------------------------------------------------------------------------------------------------
   void generateFile(const std::string& testName);

   std::vector<MeasurementValue<vfc::float32_t> > m_measurementValuesFloat;
   std::vector<MeasurementValue<bool> >           m_measurementValuesBool;
   std::vector<MeasurementValue<vfc::uint8_t> >   m_measurementValuesInt;
   std::vector<MeasurementValue<vfc::int8_t> >    m_measurementValuesSignedInt;
private:
   FILE* out;
   double       m_numberOfSignals{0};
   std::string  m_testName;
   double       m_testTime{0.0};
};

}  // end of namespace DevTest
}  // end of namespace Rb
}  // end of namespace Vmc


#endif // DC_ALCPDPDEVTEST_GENERATEMEASUREMENTFILE_HPP_INCLUDED