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
for (int i{0}; i<300; i++){
    std::string nameX = "corridorGlobalFrame_X_";
    std::string nameY = "corridorGlobalFrame_Y_";
    std::string nameOrientation = "corridorGlobalFrame_c1_";
    std::string nameCurvature = "corridorGlobalFrame_c2_";
    std::string index = std::to_string(i);
    std::string nameX_index = nameX + index;
    std::string nameY_index = nameY + index;
    std::string nameOrientation_index = nameOrientation + index;
    std::string nameCurvature_index = nameCurvature + index;
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameX_index,
        "-",
        &m_obj->corridorGlobalFrame.corridorPoints[i].PointsX.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameY_index,
        "-",
        &m_obj->corridorGlobalFrame.corridorPoints[i].PointsY.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameOrientation_index,
        "-",
        &m_obj->corridorGlobalFrame.corridorOrientation[i].value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameCurvature_index,
        "-",
        &m_obj->corridorGlobalFrame.corridorCurvature[i].value(),
        id++));
}

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "corridorGlobalFrame_Length",
        "-",
        reinterpret_cast<vfc::float32_t*>(&m_obj->corridorGlobalFrame.corridorLength),
        id++));

for (int i{0}; i<300; i++){
    std::string nameX = "corridorEgoFrame_X_";
    std::string nameY = "corridorEgoFrame_Y_";
    std::string nameOrientation = "corridorEgoFrame_c1_";
    std::string nameCurvature = "corridorEgoFrame_c2_";
    std::string index = std::to_string(i);
    std::string nameX_index = nameX + index;
    std::string nameY_index = nameY + index;
    std::string nameOrientation_index = nameOrientation + index;
    std::string nameCurvature_index = nameCurvature + index;
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameX_index,
        "-",
        &m_obj->corridor.corridorPoints[i].PointsX.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameY_index,
        "-",
        &m_obj->corridor.corridorPoints[i].PointsY.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameOrientation_index,
        "-",
        &m_obj->corridor.corridorOrientation[i].value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameCurvature_index,
        "-",
        &m_obj->corridor.corridorCurvature[i].value(),
        id++));
}

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "corridorEgoFrame_Length",
        "-",
        reinterpret_cast<vfc::float32_t*>(&m_obj->corridor.corridorLength),
        id++));

for (int i{0}; i<300; i++){
    std::string nameX = "corridorPlannerFrame_X_";
    std::string nameY = "corridorPlannerFrame_Y_";
    std::string nameOrientation = "corridorPlannerFrame_c1_";
    std::string nameCurvature = "corridorPlannerFrame_c2_";
    std::string index = std::to_string(i);
    std::string nameX_index = nameX + index;
    std::string nameY_index = nameY + index;
    std::string nameOrientation_index = nameOrientation + index;
    std::string nameCurvature_index = nameCurvature + index;
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameX_index,
        "-",
        &m_obj->corridorPlannerFrame.corridorPoints[i].PointsX.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameY_index,
        "-",
        &m_obj->corridorPlannerFrame.corridorPoints[i].PointsY.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameOrientation_index,
        "-",
        &m_obj->corridorPlannerFrame.corridorOrientation[i].value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameCurvature_index,
        "-",
        &m_obj->corridorPlannerFrame.corridorCurvature[i].value(),
        id++));
}

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "corridorPlannerFrame_Length",
        "-",
        reinterpret_cast<vfc::float32_t*>(&m_obj->corridorPlannerFrame.corridorLength),
        id++));

for (int i{0}; i<300; i++){
    std::string nameX = "trajectoryPlannerFrame_X_";
    std::string nameY = "trajectoryPlannerFrame_Y_";
    std::string index = std::to_string(i);
    std::string nameX_index = nameX + index;
    std::string nameY_index = nameY + index;
    
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameX_index,
        "-",
        &m_obj->trajectoryPlannerFrame.trajectoryPoints[i].PointsX.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameY_index,
        "-",
        &m_obj->trajectoryPlannerFrame.trajectoryPoints[i].PointsY.value(),
        id++));
}

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "trajectoryPlannerFrame_Length",
        "-",
        reinterpret_cast<vfc::float32_t*>(&m_obj->trajectoryPlannerFrame.trajectoryLength),
        id++));

for (int i{0}; i<300; i++){
    std::string nameX = "trajectoryEgoFrame_X_";
    std::string nameY = "trajectoryEgoFrame_Y_";
    std::string index = std::to_string(i);
    std::string nameX_index = nameX + index;
    std::string nameY_index = nameY + index;
    
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameX_index,
        "-",
        &m_obj->trajectoryEgoFrame.trajectoryPoints[i].PointsX.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameY_index,
        "-",
        &m_obj->trajectoryEgoFrame.trajectoryPoints[i].PointsY.value(),
        id++));
}

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "trajectoryEgoFrame_Length",
        "-",
        reinterpret_cast<vfc::float32_t*>(&m_obj->trajectoryEgoFrame.trajectoryLength),
        id++));

for (int i{0}; i<300; i++){
    std::string nameX = "trajectoryGlobalFrame_X_";
    std::string nameY = "trajectoryGlobalFrame_Y_";
    std::string index = std::to_string(i);
    std::string nameX_index = nameX + index;
    std::string nameY_index = nameY + index;
    
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameX_index,
        "-",
        &m_obj->trajectoryGlobalFrame.trajectoryPoints[i].PointsX.value(),
        id++));
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        nameY_index,
        "-",
        &m_obj->trajectoryGlobalFrame.trajectoryPoints[i].PointsY.value(),
        id++));
}

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "trajectoryGlobalFrame_Length",
        "-",
        reinterpret_cast<vfc::float32_t*>(&m_obj->trajectoryGlobalFrame.trajectoryLength),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "egoPoseGlobal_X",
        "-",
        &m_obj->egoPoseGlob.Pose2DCoordinates.PointsX.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "egoPoseGlobal_Y",
        "-",
        &m_obj->egoPoseGlob.Pose2DCoordinates.PointsY.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "egoPoseGlobal_Theta",
        "-",
        &m_obj->egoPoseGlob.Pose2DTheta.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "plannerFrame_X",
        "-",
        &m_obj->plannerFramePoseGlobalFrame.Pose2DCoordinates.PointsX.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "plannerFrame_Y",
        "-",
        &m_obj->plannerFramePoseGlobalFrame.Pose2DCoordinates.PointsY.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        "plannerFrame_Theta",
        "-",
        &m_obj->plannerFramePoseGlobalFrame.Pose2DTheta.value(),
        id++));
*/

for (int i{0}; i<7; i++){
    std::string name = "U_";
    std::string index = std::to_string(i);
    std::string name_index = name + index;
    
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        name_index,
        "-",
        &m_obj->driverModel.U[i],
        id++));
}

for (int i{0}; i<3; i++){
    std::string name = "x_";
    std::string index = std::to_string(i);
    std::string name_index = name + index;
    
    m_generateMeasurementFile->m_measurementValuesFloat.push_back(
        Dc::Emg::DevTest::MeasurementValue<vfc::float32_t>(
        name_index,
        "-",
        &m_obj->driverModel.x[i],
        id++));
}

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
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "ayLane", "-", &m_obj->ayLane.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgo", "-", &m_obj->ayEgo.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgoRelative", "-", &m_obj->ayEgoRelative.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgoRelative", "-", &m_obj->ayEgoRelative.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "ayEgoRelativeFiltered", "-", &m_obj->ayEgoRelativeFiltered.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "m_referenceLinePackage_in.c0", "-", &m_obj->m_referenceLinePackage_in.c0.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "m_referenceLinePackage_in.c1", "-", &m_obj->m_referenceLinePackage_in.c1,
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "m_referenceLinePackage_in.c2", "-", &m_obj->m_referenceLinePackage_in.c2.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "m_pdpLaneChangeDetection.m_laneChangeLateralDisplacement", "-", &m_obj->m_pdpLaneChangeDetection.m_laneChangeLateralDisplacement.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "m_pdpLaneChangeDetection.m_laneChangeDuration", "-", &m_obj->m_pdpLaneChangeDetection.m_laneChangeDuration.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
        "m_pdpLaneChangeDetection.m_laneChangeLateralSpeed", "-", &m_obj->m_pdpLaneChangeDetection.m_laneChangeLateralSpeed.value(),
        id++));

m_generateMeasurementFile->m_measurementValuesFloat.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::float32_t>(
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
	Rb::Vmc::DevTest::MeasurementValue<vfc::uint8_t>(
		"m_pdpLaneChangeDetection.m_laneChangeState",
		"",
		reinterpret_cast<vfc::uint8_t*>(&m_obj->m_pdpLaneChangeDetection.m_laneChangeState),
		id++,
		"{0=\"NoLaneChange\",1=\"LaneChangeUnverified\",2=\"LaneChangeVerified\"}"));

m_generateMeasurementFile->m_measurementValuesInt.push_back(
	Rb::Vmc::DevTest::MeasurementValue<vfc::uint8_t>(
		"m_pdpLaneChangeDetection.m_laneChangeDirection",
		"",
		reinterpret_cast<vfc::uint8_t*>(&m_obj->m_pdpLaneChangeDetection.m_laneChangeDirection),
		id++,
		"{0=\"pdp_LaneChangeUnknown\",1=\"pdp_LaneChangeLeft\",2=\"pdp_LaneChangeRight\"}"));*/
