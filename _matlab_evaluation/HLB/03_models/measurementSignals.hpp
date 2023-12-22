/*
 * Measurement signal list for d97
 * Original functionality by Schmid Steffen (XC-DX/EYM5)
 * Add desired measurement variables in ordered, type-specific blocks:
 * 1. floats
 * 2. enums
 * 3. bools
 * Use following type-specific templates:
*/

/*
  //------------------------------floats-----------------------------------
m_generateMeasurementFile->m_measurementValuesFloat.push_back(
Dc::Mcam::Lat::DevTest::MeasurementValue<vfc::float32_t>(
  "m_sacControlCascade.m_steerAngleSpeedControl.m_finalTorque.m_steerTqTarTBTLimited",
  "-",
  &m_obj->m_sacControlCascade.m_steerAngleSpeedControl.m_finalTorque.m_steerTqTarTBTLimited.value(),
  id++));

//------------------------------bools-----------------------------------
std::string enumspecBit = "{0=\"false\", 1=\"true\"}";
m_generateMeasurementFile->m_measurementValuesBool.push_back(
   Dc::Mcam::Lat::DevTest::MeasurementValue<bool>(
      "m_sacControlCascade.m_steerAngleSpeedControl.m_PidControl.m_isCtrlLimited",
      "",
      &m_obj->m_sacControlCascade.m_steerAngleSpeedControl.m_PidControl.m_isCtrlLimited,
      id++,
      enumspecBit));
//------------------------------enums-----------------------------------
m_generateMeasurementFile->m_measurementValuesInt.push_back(
   Dc::Mcam::Lat::DevTest::MeasurementValue<vfc::uint8_t>(
      "EpsRequest",
      "",
      reinterpret_cast<vfc::uint8_t*>(&m_obj->m_McamLatSlowEpsRequestMeasLowLvl),
      id++,
      "{0=\"NoReq\", 1=\"Lvl2_TqIf\", 2=\"Lvl2_AgIf\", 3=\"Lvl3\"}"));
*/

// Start signal list

//------------------------------1. floats-----------------------------------

/*m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "ayLane", "-", &m_obj->ayLane.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgo", "-", &m_obj->ayEgo.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgoRelative", "-", &m_obj->ayEgoRelative.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgoRelative", "-", &m_obj->ayEgoRelative.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgoRelativeFiltered", "-", &m_obj->ayEgoRelativeFiltered.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "m_referenceLinePackage_in.c0", "-", &m_obj->m_referenceLinePackage_in.c0.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "m_referenceLinePackage_in.c1", "-", &m_obj->m_referenceLinePackage_in.c1,
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "m_referenceLinePackage_in.c2", "-", &m_obj->m_referenceLinePackage_in.c2.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "m_pdpLaneChangeDetection.m_laneChangeLateralDisplacement", "-", &m_obj->m_pdpLaneChangeDetection.m_laneChangeLateralDisplacement.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "m_pdpLaneChangeDetection.m_laneChangeDuration", "-", &m_obj->m_pdpLaneChangeDetection.m_laneChangeDuration.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "m_pdpLaneChangeDetection.m_laneChangeLateralSpeed", "-", &m_obj->m_pdpLaneChangeDetection.m_laneChangeLateralSpeed.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "m_pdpLaneChangeDetection.m_laneChangeFinishDebounceTimer", "-", &m_obj->m_pdpLaneChangeDetection.m_laneChangeFinishDebounceTimer.value(),
        id++));

//------------------------------2. bools-----------------------------------
//std::string enumspecBit = "{0=\"false\", 1=\"true\"}";
//m_generateMeasurementFile->m_measurementValuesBool.push_back(
//	Dc::Mcam::Lat::DevTest::MeasurementValue<bool>(
//		"m_vehicleSteerCharacteristics.m_isEstimatorActive",
//		"",
//		reinterpret_cast<bool*>(&m_obj->m_vehicleSteerCharacteristics.m_isEstimatorActive),
//		id++,
//		enumspecBit));

//------------------------------3. enums-----------------------------------
m_generateMeasurementFile->m_measurementValuesInt.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::uint8_t>(
		"m_pdpLaneChangeDetection.m_laneChangeState",
		"",
		reinterpret_cast<vfc::uint8_t*>(&m_obj->m_pdpLaneChangeDetection.m_laneChangeState),
		id++,
		"{0=\"NoLaneChange\",1=\"LaneChangeUnverified\",2=\"LaneChangeVerified\"}"));

m_generateMeasurementFile->m_measurementValuesInt.push_back(
	Dc::Emg::DevTest::MeasurementValue<vfc::uint8_t>(
		"m_pdpLaneChangeDetection.m_laneChangeDirection",
		"",
		reinterpret_cast<vfc::uint8_t*>(&m_obj->m_pdpLaneChangeDetection.m_laneChangeDirection),
		id++,
		"{0=\"pdp_LaneChangeUnknown\",1=\"pdp_LaneChangeLeft\",2=\"pdp_LaneChangeRight\"}"));*/
