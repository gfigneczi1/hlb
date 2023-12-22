# Python program to read
# json file


import json
import os
import pandas
import pandas as pd


class configurator():
    def __init__(self, inputs, wd):
        self.wd = wd
        self.configpool = os.path.join(self.wd, '_configuration\configpool.json')
        self.input = inputs

    def filesearcher(self):
        files = os.listdir(self.input['measurement_folder'])
        extensions = []
        numberOfFormat = 0
        possible_ext = 'na'
        for file in files:
            filename, file_extension = os.path.splitext(file)
            extensions.append(file_extension)
        extensions_df = pd.DataFrame(extensions)
        for possible_ext in self.input['possible file formats']:
            numberOfFormatTemp = len(extensions_df[extensions_df[0] == ('.'+possible_ext)])
            if numberOfFormatTemp > numberOfFormat:
                numberOfFormat = numberOfFormatTemp
                extension = possible_ext
        outfiles = []
        for file in files:
            if extension in file:
                file = os.path.join(self.input["measurement_folder"], file)
                outfiles.append(file)

        return outfiles, extension

    def configmapper(self, configpool):
        currentconfig = {}
        files, file_extension = self.filesearcher()

        currentconfig['extension'] = file_extension
        currentconfig['measurement_files'] = files
        signal_profiles = configpool['signal_profiles']
        signal_profiles = signal_profiles[0]
        for signal_profile, values in signal_profiles.items():
            if self.input['vehicle_config'] in signal_profile:
                if (file_extension in values):
                    currentconfig['signal_list'] = values[currentconfig['extension']]
                else:
                    currentconfig['signal_list'] = ['na']
                break
        if not currentconfig['signal_list']:
            currentconfig['signal_list'] = ['nan']
        for profile in configpool['segmentation_profiles']:
            if profile == self.input['segmentation_profile']:
                currentconfig['segmentation_profile'] = profile
        if not currentconfig['segmentation_profile']:
            currentconfig['segmentation_profile'] = 'na'

        for profile in configpool['evaluation_profiles']:
            if profile == self.input['evaluation_profile']:
                currentconfig['evaluation_profile'] = profile
        if not currentconfig['evaluation_profile']:
            currentconfig['evaluation_profile'] = 'na'

        currentconfig['vehicle_profile'] = self.input['vehicle_config']
        currentconfig['training_step'] = self.input['training_step']
        currentconfig['gentle_threshold'] = self.input['gentle_threshold']
        currentconfig['strict_threshold'] = self.input['strict_threshold']
        currentconfig['classification_model'] = self.input['classification_model']
        currentconfig['regression_model'] = self.input['regression_model']

        return currentconfig

    def execute(self):
        # Opening JSON config pool file
        fconfigpool = open(self.configpool)
        configpool = json.load(fconfigpool)

        # Closing file
        fconfigpool.close()

        #Generating current config
        currentconfig = self.configmapper(configpool)

        return currentconfig