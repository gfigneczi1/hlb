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

#include "alcpdpdevtest_generateMeasurementFile.hpp"

namespace Dc
{
namespace Emg
{
namespace DevTest
{

void GenerateMeasurementFile::writeDataD97(const vfc::CSI::si_second_f32_t& cycleTime)
{
	// data section <m>
	// ------------------------------------------------------
	// Number of signal values
	fwrite(&m_numberOfSignals, sizeof(double), 1, out);
	// Output data
	fwrite(&m_testTime, sizeof(double), 1, out);
	m_testTime += cycleTime.value();
   // Write Float data
	for (std::vector<MeasurementValue<vfc::float32_t> >::iterator v = std::begin(m_measurementValuesFloat);
		v != std::end(m_measurementValuesFloat);
		++v)
	{
	  // Convert to double and write binary data
	  double val = *v->valueAddress;
	  fwrite(&val, sizeof(double), 1, out);
	}
   // Write Bool data
	for (std::vector<MeasurementValue<bool> >::iterator v = std::begin(m_measurementValuesBool);
		v != std::end(m_measurementValuesBool);
		++v)
	{
	  // Convert to double and write binary data
	  double val = *v->valueAddress;
	  fwrite(&val, sizeof(double), 1, out);
	}
   // Write Enum(Int) data
   for (std::vector<MeasurementValue<vfc::uint8_t> >::iterator v = std::begin(m_measurementValuesInt);
        v != std::end(m_measurementValuesInt);
        ++v)
   {
      // Convert to double and write binary data
      double val = *v->valueAddress;
	   fwrite(&val, sizeof(double), 1, out);
   }
    // Write Signed Int data
   for (std::vector<MeasurementValue<vfc::int8_t> >::iterator v = std::begin(m_measurementValuesSignedInt);
        v != std::end(m_measurementValuesSignedInt);
        ++v)
   {
      // Convert to double and write binary data
      double val = *v->valueAddress;
	   fwrite(&val, sizeof(double), 1, out);
   }
}

void GenerateMeasurementFile::closeFile()
{
   fclose(out);
}

void GenerateMeasurementFile::generateFile(const std::string& testName)
{   
    char stringRes[5000000]; // buffer string array for signal headers, this needs to be "large enough"
	//out = fopen("test.d97", "wb");
	out = fopen((testName +std::string(".d97")).c_str(), "wb");
	m_numberOfSignals = m_measurementValuesFloat.size() + m_measurementValuesBool.size() + m_measurementValuesInt.size() + m_measurementValuesSignedInt.size() + 1.0;	
	
	// Write header
	sprintf(stringRes, "[CONTROL]\nFLOAT=INTEL387\nPRECISION=DOUBLE\nDATAMATRIX=FRAME_BY_*\nTMSV=Time\n");
	sprintf(stringRes+strlen(stringRes), "VDIM=%d\n\n", uint32_t(m_numberOfSignals));

	// Define signals
	// time vector
   sprintf(stringRes+strlen(stringRes), "[SIGNAL0]\nNAME=Time\nDPOS=0\nDESCR=TIME BASE\nUNIT=s\n\n");
	// All other signals
	// Float signals
	for (std::vector<MeasurementValue<vfc::float32_t> >::iterator v = std::begin(m_measurementValuesFloat);
		v != std::end(m_measurementValuesFloat);
		++v)
	{
		sprintf(stringRes+strlen(stringRes), "[SIGNAL%d]\nNAME=%s\nDPOS=%d\nDESCR=TIME BASE\nUNIT=%s\n\n", v->id, v->name.c_str(), v->id, v->unit.c_str());
	}
	// Bool signals
	for (std::vector<MeasurementValue<bool> >::iterator v = std::begin(m_measurementValuesBool);
		v != std::end(m_measurementValuesBool);
		++v)
	{
		sprintf(stringRes+strlen(stringRes), "[SIGNAL%d]\nNAME=%s\nDPOS=%d\nDESCR=TIME BASE\nUNIT=%s\n\n", v->id, v->name.c_str(), v->id, v->unit.c_str());
	}
   // Enum(Int) signals
   for (std::vector<MeasurementValue<vfc::uint8_t> >::iterator v = std::begin(m_measurementValuesInt);
        v != std::end(m_measurementValuesInt);
        ++v)
   {
      sprintf(stringRes+strlen(stringRes), "[SIGNAL%d]\nNAME=%s\nDPOS=%d\nDESCR=TIME BASE\nUNIT=%senum=%s\n\n", v->id, v->name.c_str(), v->id, v->unit.c_str(), v->enumspec.c_str());
   }
   
    // Enum(Signed Int) signals
   for (std::vector<MeasurementValue<vfc::int8_t> >::iterator v = std::begin(m_measurementValuesSignedInt);
        v != std::end(m_measurementValuesSignedInt);
        ++v)
   {
      sprintf(stringRes+strlen(stringRes), "[SIGNAL%d]\nNAME=%s\nDPOS=%d\nDESCR=TIME BASE\nUNIT=%senum=%s\n\n", v->id, v->name.c_str(), v->id, v->unit.c_str(), v->enumspec.c_str());
   }
   // Write data
	sprintf(stringRes+strlen(stringRes), "[DATA]\n");
	fwrite(&stringRes, sizeof(char), strlen(stringRes), out);
	
	// Number of rows in the data section <n> -> Frames_By setting: Number of samples
	// Number of samples is not important with Frames_By setting --> set to 0.0
	double samples = 0;
	fwrite(&samples, sizeof(samples), 1, out);


	//fclose(out);
}


}  // end of namespace DevTest
}  // end of namespace Emg
}  // end of namespace Dc
