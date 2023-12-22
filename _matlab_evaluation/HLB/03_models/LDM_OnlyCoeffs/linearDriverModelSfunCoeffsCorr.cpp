#define _ALLOW_KEYWORD_MACROS
#define private public
#define protected public
#include "emg_linearDriverModel.hpp"
#include "./alcpdpdevtest_generateMeasurementFile.hpp"

#undef private
#undef protected
#undef _ALLOW_KEYWORD_MACROS

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME emg_linearDriverModelCoeffsCorr

// user parameters
#define NINP 11          // number of Simulink-level inputs (without parameter port)
#define NOUP 2          // number of Simulink-level outputs
#define SAMPLETIME 0.02  // component (sfun) sample time
#define DATALOGGING true // flag for data logging into file - faster simtime without

// Need to include simstruc.h for the definition of the SimStruct and
// its associated macro definitions.
#include "simstruc.h"

// #define IS_PARAM_DOUBLE(pVal)                                                                              \
//    (mxIsNumeric(pVal) && !mxIsLogical(pVal) && !mxIsEmpty(pVal) && !mxIsSparse(pVal) && !mxIsComplex(pVal) \
//     && mxIsDouble(pVal))

// Function: mdlInitializeSizes ===============================================
// Abstract:
//    The sizes information is used by Simulink to determine the S-function
//    block's characteristics (number of inputs, outputs, states, etc.).
static void mdlInitializeSizes(SimStruct *S)
{
   // No expected parameters
   ssSetNumSFcnParams(S, 0);

   // Parameter mismatch will be reported by Simulink
   if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
   {
      return;
   }

   // Specify I/O
   if (!ssSetNumInputPorts(S, NINP))
      return;
   for (int ii = 0; ii < NINP; ii++)
   { // set scalar inputs
      if (ii < 4){
          ssSetInputPortWidth(S, ii, 1);
      }
      else if (ii == 9){
          ssSetInputPortWidth(S, ii, 21);
      }
      else if (ii == 10){
          ssSetInputPortWidth(S, ii, 3);
      }
      else{
          ssSetInputPortWidth(S, ii, 1);
      }
      ssSetInputPortDirectFeedThrough(S, ii, 1);
   }
   if (!ssSetNumOutputPorts(S, NOUP))
      return;
   for (int ii = 0; ii < NOUP; ii++)
   { // assumption: only scalar outputs
        ssSetOutputPortWidth(S, ii, 30);  
   }

   ssSetNumSampleTimes(S, 1);

   // Reserve place for C++ objects
   ssSetNumPWork(S, 3);

   // ssSetOperatingPointCompliance(S, USE_CUSTOM_OPERATING_POINT);

   // ssSetOptions(
   //    S,
   //    SS_OPTION_WORKS_WITH_CODE_REUSE | SS_OPTION_EXCEPTION_FREE_CODE
   //       | SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME);
}

// Function: mdlInitializeSampleTimes =========================================
// Abstract:
//   This function is used to specify the sample time(s) for your
//   S-function. You must register the same number of sample times as
//   specified in ssSetNumSampleTimes.
static void mdlInitializeSampleTimes(SimStruct *S)
{
   ssSetSampleTime(S, 0, SAMPLETIME);
   ssSetOffsetTime(S, 0, 0.0);
   ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

// Function: mdlStart =======================================================
// Abstract:
//   This function is called once at start of model execution. If you
//   have states that should be initialized once, this is the place
//   to do it.
#define MDL_START
static void mdlStart(SimStruct *S)
{
   // Store new C++ object in the pointers vector
   Dc::Emg::LinearDriverModel* m_obj = new Dc::Emg::LinearDriverModel ();
   ssGetPWork(S)[0] = m_obj;

   // setup measurement config
   Dc::Emg::DevTest::GenerateMeasurementFile *m_generateMeasurementFile =
       new Dc::Emg::DevTest::GenerateMeasurementFile();

   int id = 1;
#include "measurementSignals.hpp"
   // create measurement file and store handle
   m_generateMeasurementFile->generateFile("driverModel");
   ssGetPWork(S)[2] = m_generateMeasurementFile;
}

// Function: mdlOutputs =======================================================
// Abstract:
//   In this function, you compute the outputs of your S-function
//   block.
static void mdlOutputs(SimStruct *S, int_T tid)
{
   // Retrieve C++ object from the pointers vector
   Dc::Emg::LinearDriverModel *m_obj = reinterpret_cast<Dc::Emg::LinearDriverModel *>(ssGetPWork(S)[0]);
   
   Dc::Emg::DevTest::GenerateMeasurementFile *m_generateMeasurementFile =
       reinterpret_cast<Dc::Emg::DevTest::GenerateMeasurementFile *>(ssGetPWork(S)[2]);
   
   // Get data addresses of I/O
   int inputIndex = 0;
   InputRealPtrsType c0 = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType c1 = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType c2 = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType c3 = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType egoPoseX = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType egoPoseY = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType egoPoseTheta = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType LAD = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType replanCycle = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType P = ssGetInputPortRealSignalPtrs(S, inputIndex++);
   InputRealPtrsType P_nodePointDistances = ssGetInputPortRealSignalPtrs(S, inputIndex++);

   // get output pointers
   int outputIndex = 0;
   real_T *trajectoryCoeffsThreeSegments = ssGetOutputPortRealSignal(S, outputIndex++);
   real_T *nodePointsEgo = ssGetOutputPortRealSignal(S, outputIndex++);
   
   vfc::float32_t P_in[21]{0.0f};
   vfc::float32_t P_nodePointDistances_in[3]{0.0f};
   
   Dc::Emg::CorridorInfoCoefficients corridorCoefficients;
   Dc::Emg::Pose2D egoPoseGlobal;
   Dc::Emg::LDMParamIn parameters;
   
   for (int i{0}; i < 21; i++){
       parameters.P[i] = static_cast<vfc::float32_t>(*P[i]);
   }
   
   for (int i{0}; i < 3; i++){
       parameters.P_nodePointDistances[i] = static_cast<vfc::float32_t>(*P_nodePointDistances[i]);
   }   
   
   corridorCoefficients.c0 = static_cast<vfc::float32_t>(*c0[0]);
   corridorCoefficients.c1 = static_cast<vfc::float32_t>(*c1[0]);
   corridorCoefficients.c2 = static_cast<vfc::float32_t>(*c2[0]);
   corridorCoefficients.c3 = static_cast<vfc::float32_t>(*c3[0]);
   
   egoPoseGlobal.Pose2DCoordinates.PointsX = static_cast<vfc::CSI::si_metre_f32_t>(*egoPoseX[0]);
   egoPoseGlobal.Pose2DCoordinates.PointsY = static_cast<vfc::CSI::si_metre_f32_t>(*egoPoseY[0]);
   egoPoseGlobal.Pose2DTheta = static_cast<vfc::CSI::si_radian_f32_t>(*egoPoseTheta[0]);
   
   parameters.lookAheadDistance = static_cast<vfc::CSI::si_metre_f32_t>(*LAD[0]);
   parameters.replanCycle = static_cast<vfc::uint8_t>(*replanCycle[0]);
   
   Dc::Emg::PolynomialCoeffs5TwoSegments y = m_obj->runCoeffsLite(
           corridorCoefficients,
           egoPoseGlobal,
           parameters);
   
   for (int i{0}; i < 2; i++)
   {
       trajectoryCoeffsThreeSegments[0+(i*8)] = y.segmentCoeffs[i].c0;
       trajectoryCoeffsThreeSegments[1+(i*8)] = y.segmentCoeffs[i].c1;
       trajectoryCoeffsThreeSegments[2+(i*8)] = y.segmentCoeffs[i].c2;
       trajectoryCoeffsThreeSegments[3+(i*8)] = y.segmentCoeffs[i].c3;
       trajectoryCoeffsThreeSegments[4+(i*8)] = y.segmentCoeffs[i].c4;
       trajectoryCoeffsThreeSegments[5+(i*8)] = y.segmentCoeffs[i].c5;
       trajectoryCoeffsThreeSegments[6+(i*8)] = y.segmentCoeffs[i].length.value();
       trajectoryCoeffsThreeSegments[7+(i*8)] = y.segmentCoeffs[i].breakPoint.value();
   }
   
   for (int i{0}; i<4; i++)
   {
       nodePointsEgo[0+(i*2)] = m_obj->nodePointsEgoFrame.nodePointsCoordinates[i].PointsX.value();
       nodePointsEgo[1+(i*2)] = m_obj->nodePointsEgoFrame.nodePointsCoordinates[i].PointsY.value();
   }

   if (DATALOGGING)
   {
      m_generateMeasurementFile->writeDataD97(vfc::CSI::si_second_f32_t {0.02f});
   }
}
// #ifdef MATLAB_MEX_FILE
// /* Define to indicate that this S-Function has the mdlG[S]etOperatingPoint methods */
// #define MDL_OPERATING_POINT
//
// /* Function: mdlGetOperatingPoint ==================================================
//  * Abstract:
//  *    Save the operating point of this block and return it to Simulink
//  */
// static mxArray* mdlGetOperatingPoint(SimStruct* S)
// {
//     DoubleAdder *da = static_cast<DoubleAdder*>(ssGetPWork(S)[0]);
//     return mxCreateDoubleScalar(da->GetPeak());
// }
// /* Function: mdlSetOperatingPoint =================================================
//  * Abstract:
//  *   Restore the operating point of this block based on the provided data (ma)
//  *   The data was saved by mdlGetOperatingPoint
//  */
// static void mdlSetOperatingPoint(SimStruct* S, const mxArray* ma)
// {
//     // Retrieve C++ object from the pointers vector
//     DoubleAdder *da = static_cast<DoubleAdder*>(ssGetPWork(S)[0]);
//     da->SetPeak(mxGetPr(ma)[0]);
// }
// #endif // MATLAB_MEX_FILE

// Function: mdlTerminate =====================================================
// Abstract:
//   In this function, you should perform any actions that are necessary
//   at the termination of a simulation.  For example, if memory was
//   allocated in mdlStart, this is the place to free it.
static void mdlTerminate(SimStruct *S)
{
   // Retrieve and destroy C++ object
   Dc::Emg::LinearDriverModel  *m_obj = reinterpret_cast<Dc::Emg::LinearDriverModel  *>(ssGetPWork(S)[0]);
   delete m_obj;
   Dc::Emg::DevTest::GenerateMeasurementFile *m_generateMeasurementFile =
       reinterpret_cast<Dc::Emg::DevTest::GenerateMeasurementFile *>(ssGetPWork(S)[2]);
   m_generateMeasurementFile->closeFile();
   delete m_generateMeasurementFile;
}

// Required S-function trailer
#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */
#include "simulink.c"  /* MEX-file interface mechanism */
#else
#include "cg_sfun.h" /* Code generation registration function */
#endif
