# Preprocess components:
# 1. configuration
# 1.1 Input read in
# 1.2 configuration pool read in
# 2. measurement conversion
# 2.1 measurement file read-in
# 2.2 interface file generation

import os.path
import shutil
from os import path
import json

from _configuration.configuration import configurator
from _conversion_handler.conversion_handler import Convert
from gooey import GooeyParser, Gooey
import os
import argparse

class preprocess():
    def __init__(self):
        self.rootfolder = os.path.dirname(os.path.abspath(__file__))

    def tempfolderhandler(self):
        if path.isdir(path.join(self.rootfolder, '_temp')):
            shutil.rmtree(path.join(self.rootfolder, '_temp'))

        os.mkdir(path.join(self.rootfolder, '_temp'))

    @Gooey(program_name="Input Handler", default_size=(600, 600), tabbed_groups=True)
    def execute(self):
        # STEP 0: temp folder handling
        self.tempfolderhandler()
        # STEP 1: Input reading
        with open(path.join(self.rootfolder, '_configuration', 'configpool.json')) as f:
            configpool = json.load(f)
        vehicle_config = [x.split(sep="_") for x in tuple(configpool["signal_profiles"][0].keys())]
        vehicle_config = ["_".join(s[2:]) for s in vehicle_config]
        if path.exists(path.join(self.rootfolder, '_configuration', 'previous_inputs.json')):
            with open(path.join(self.rootfolder, '_configuration', 'previous_inputs.json')) as f:
                previous_inputs = json.load(f)
        else:
            previous_inputs = {
                "measurement_folder": None,
                "evaluation_profile": None,
                "vehicle_config": None,
                "segmentation_profile": None,
                "delete_temp_files": None,
                "gentle_threshold": None,
                "strict_threshold": None
            }
            with open(path.join(self.rootfolder, '_configuration', 'previous_inputs.json'), 'w', encoding='utf-8') as f:
                json.dump(previous_inputs, f, ensure_ascii=False, indent=5)

        parser = GooeyParser(description="Parsing inputs for configurator, only .mf4, .bag and"
                                         " .mat file extensions are supported")
        main_group = parser.add_argument_group("Main Options","Select Main Options")
        main_group.add_argument("-i", "--MeasurementFolder", required=True, default=previous_inputs["measurement_folder"],
                            help="Path to the measurement folder", widget='DirChooser')
        main_group.add_argument("-e", "--EvaluationProfile", required=True, default=previous_inputs["evaluation_profile"],
                            help="Select the evaluation profile", choices=configpool["evaluation_profiles"],
                            widget='Dropdown')
        main_group.add_argument("-v", "--VehicleConfig", required=True, default=previous_inputs["vehicle_config"],
                            help="Select the vehicle config", choices=vehicle_config,
                            widget='Dropdown')
        main_group.add_argument("-s", "--SegmentationProfile", required=True, default=previous_inputs["segmentation_profile"],
                            help="Select the segmentation profile", choices=configpool["segmentation_profiles"],
                            widget='Dropdown')
        main_group.add_argument("-d", "--DeleteTempFiles", default=previous_inputs["delete_temp_files"],
                            help="If True, deletes the content of the _temp folder after new run (default: False)",
                            action='store_true')
        py_group = parser.add_argument_group('PythonOptions', "Select Python Options. Only apply if an Evaluation profile with the word 'Python' is selected.")
        eval_group = py_group.add_argument_group("Evaluation Options", "Only apply if Evaluation profile includes Python and is not train.")
        train_group = py_group.add_argument_group("Training Options", "Only apply if 'Python_train' is selected as the Evaluation profile.")
        py_group.add_argument("-c", "--ClassificationModel", default="Classification_ResNet-20_100",
                            help="Select which classification model to use.",
                            choices=["Classification_ResNet-20_100", "class_model1", "class_model2", "class_model3"],
                            widget='Dropdown')
        py_group.add_argument("-r", "--RegressionModel", default="Simple_Regression_model",
                            help="Select which regression model to use.", choices=["Simple_Regression_model", "reg_model1", "reg_model2", "reg_model3"],
                            widget='Dropdown')
        eval_group.add_argument("-g", "--GentleThreshold", default=previous_inputs["gentle_threshold"],
                            help="Allowed gentle Margin of Error in seconds.",
                            widget='DecimalField')
        eval_group.add_argument("-t", "--StrictThreshold", default=previous_inputs["strict_threshold"],
                            help="Allowed strict Margin of Error in seconds.",
                            widget='DecimalField')
        train_group.add_argument("-y", "--TrainingStep", default="preprocessing", help="Select which training step will happen.",
                            choices=["preprocessing", "train_classification", "train_regression", "cleanup"], widget='Dropdown')
        args = parser.parse_args()
        inputs = {'measurement_folder': args.MeasurementFolder,
                  'evaluation_profile': args.EvaluationProfile,
                  'vehicle_config': args.VehicleConfig,
                  'delete_temp_files': str(args.DeleteTempFiles).lower(),
                  'possible file formats': ['mf4', 'bag', 'mat'],
                  'segmentation_profile': args.SegmentationProfile,
                  'training_step': args.TrainingStep,
                  'gentle_threshold': args.GentleThreshold,
                  'strict_threshold': args.StrictThreshold,
                  'classification_model': args.ClassificationModel,
                  'regression_model': args.RegressionModel
                  }

        # inputs = {'measurement_folder': r'C:\database\TSS4_PBP_ALC_2022_07\base_learning\toEval',
        #           'evaluation_profile': "Python_evaluation_online",
        #             'vehicle_config': "signal_profile_VW_Golf_Bp_2211_ES910_CAN_only",
        #             'possible file formats': ['mf4', 'bag', 'mat'],
        #             'segmentation_profile': "Resimulation",
        #             'training_step': "preprocessing",
        #             'gentle_threshold': 1.5,
        #             'strict_threshold': 0.5,
        #             'classification_model': "Classification_ResNet-20_100",
        #             'regression_model': "Simple_Regression_model"
        #           }

        previous_inputs["measurement_folder"] = args.MeasurementFolder
        previous_inputs["evaluation_profile"] = args.EvaluationProfile
        previous_inputs["vehicle_config"] = args.VehicleConfig
        previous_inputs["segmentation_profile"] = args.SegmentationProfile
        previous_inputs["delete_temp_files"] = args.DeleteTempFiles
        previous_inputs["gentle_threshold"] = args.GentleThreshold
        previous_inputs["strict_threshold"] = args.StrictThreshold
        with open(path.join(self.rootfolder, '_configuration', 'previous_inputs.json'), 'w', encoding='utf-8') as f:
            json.dump(previous_inputs, f, ensure_ascii=False, indent=5)
        # STEP 2: Configuration
        configurator_engine = configurator(inputs, self.rootfolder)
        currentconfig = configurator_engine.execute()

        with open(path.join(self.rootfolder, '_temp', 'config.json'), 'w', encoding='utf-8') as f:
            json.dump(currentconfig, f, ensure_ascii=False, indent=4)

        # STEP 3: Convert MDF Files to .mat format
        Convert(currentconfig["extension"], currentconfig["measurement_files"], self.rootfolder,
                currentconfig["signal_list"], currentconfig["vehicle_profile"]).convert_handler()

        return currentconfig

