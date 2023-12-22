# Acceleration profiler: an ACC's acceleration personalization prototype
The personalization of the ACC's acceleration project is a personalization functionality prototype to evaluate the possibility of adapting the ACC's acceleration to the driver's style. The main goals are to identify the driving style of a driver and adapt the ACC to the style.

This project aims to analyze vehicle data to deduce driver's behavior during accelerations.

Contents
- [General information](#general-information)
- [Requirements](#requirements)
- [Python environment setup](#python-environment-setup)
    - [Conda](#conda)
- [Profiler usage](#profiler-usage)
    - [Help](#help)
    - [Configuration](#configuration)
    - [Graphical interface](#graphical-interface)
    - [Command line interface](#command-line-interface)
- [Analysis prototypes](#analysis-prototypes)
    - [Models' parameters](#models-parameters)
    - [Weighted formula tuning](#weighted-formula-tunning)
    - [Speedup curvature metric](#speedup-curvature-metric)
- [Speedup and lane change visualization](#speedup-and-lane-change-visualization)
- [Testing](#testing)
    - [Time window extraction](#time-window-extraction)
## General information
- Source code: [sourcecode01.de.bosch.com/projects/VMCBPTEAM/repos/kdp_hlb_evalframework](https://sourcecode01.de.bosch.com/projects/VMCBPTEAM/repos/kdp_hlb_evalframework)
- Documentation: [inside-docupedia.bosch.com/confluence/display/ECL/%5BLiottard+Julien+and+Schoenenberger+Samuel+Naoki%5D+Acceleration+Learning+for+ACC+personalization](https://inside-docupedia.bosch.com/confluence/display/ECL/%5BLiottard+Julien+and+Schoenenberger+Samuel+Naoki%5D+Acceleration+Learning+for+ACC+personalization)
- Curated vehicle data per trip: [\\\abtvdfs2.de.bosch.com\ismdfs\ilm\abt\AVT_Datasets\CAN_Exports\\](\\\abtvdfs2.de.bosch.com\ismdfs\ilm\abt\AVT_Datasets\CAN_Exports\\) or [https://inside-docupedia.bosch.com/confluence/display/AVT/Curated+timeseries+CAN+data](https://inside-docupedia.bosch.com/confluence/display/AVT/Curated+timeseries+CAN+data)
- Subject mapping to trip records: [https://sourcecode.socialcoding.bosch.com/projects/TEAM_ADAS/repos/fdas/browse/extractlabelfeatures/processing/scripts/ArchiveCheck_update.csv](https://sourcecode.socialcoding.bosch.com/projects/TEAM_ADAS/repos/fdas/browse/extractlabelfeatures/processing/scripts/ArchiveCheck_update.csv) or [\\\abtvdfs2.de.bosch.com\ismdfs\ilm\abt\AVT_Datasets\trip_keys](\\\abtvdfs2.de.bosch.com\ismdfs\ilm\abt\AVT_Datasets\trip_keys)
## Requirements
- Python 3.7
- Python 3.10 for map data extraction (optional)
- Pip 22.3.1
- Anaconda 4.11
- Git 2.40
## Python environment setup
### Conda
To create a conda environment and install the Python packages, run:
```bash
conda create --name acceleration_profiler_env python=3.7 -y
conda activate acceleration_profiler_env
pip install -r requirements.txt
```
To exit the newly created conda environment, run:
```bash
conda deactivate
```
## Profiler usage
The software architecture is described in the documentation (see [General information](#general-information)). It can extract time windows of the data and estimate parameters of models.
### Help
```bash
python acceleration_profiler_main.py -h
```
### Configuration
You can configure the paths to the different files and folders used by the profiler by editing the TOML configuration file. The default configuration file used by the program is `_python_evaluation\configuration\default_paths_configuration.toml`
### Graphical interface
For information, you can read the following Docupedia page: [https://inside-docupedia.bosch.com/confluence/display/ECL/The+Acceleration+Profile-inator](https://inside-docupedia.bosch.com/confluence/display/ECL/The+Acceleration+Profile-inator). You can enter the parameters using a graphical interface by running:
```bash
python acceleration_profiler_main.py --gui
```
### Command line interface
You can also enter the parameters with the command line interface.
For example, you can extract time window from a single trip data file, by running:
```bash
python acceleration_profiler_main.py --step preprocessing --InputMode File --MappingFile C:\path\to\subject_trip_mapping_file.csv --InputFile C:\path\to\trip_data_file.csv
```
For example, you can estimate CA PT2 model's parameters from a single file, by running:
```bash
python acceleration_profiler_main.py --step modeling --InputMode File --MappingFile C:\path\to\subject_trip_mapping_file.csv --InputFile C:\path\to\trip_data_file.csv --Estimation_Step python_estimation -model CA_PT2 --MATLAB_command parameter_estimation --Starting_Speedup_Index 0 --Ending_Speedup_Index 1 --Save_postfix version1
```
## Analysis prototypes
### Models' parameters
You can view the results of the analysis and clustering of the models (PT1, PT2 and CA-PT2) in the following Jupyter Notebook located in the `acceleration_learning\analysis` folder:
```
PT1_parameter_analysis.ipynb
PT2_parameter_analysis.ipynb
CA_PT2_parameter_analysis.ipynb
```
### Weighted formula tunning
You can run an example of use of the weights formula tuner module:
```bash
python weights_tuning_main.py
```
### Speedup curvature metric
You can run an example of use of the speedup curvature module:
```bash
python curvature_analysis_main.py
```
## Troubleshooting
### Speedup and lane change visualization
A tool was developed to troubleshoot the speedup extraction function. This tool displays velocity data of a trip in time for speedup and the lane marking position for lane change.
For example to display the velocity of a trip data by having the extracted speedup highlighted, run:
```bash
python visualize_situation_windows.py --directory C:\path\to\trip_data_directory --scenario speed_up --filename trip_data_file.csv
```
### Video feedback from trip data
A Python script to play the video feedback displaying some signals' values can be used with the visualization to debug. The script is available at [https://sourcecode.socialcoding.bosch.com/projects/TEAM_ADAS/repos/fdal/browse/ExplorationTools/VideoPlayer](https://sourcecode.socialcoding.bosch.com/projects/TEAM_ADAS/repos/fdal/browse/ExplorationTools/VideoPlayer), for example `ENA1_AVT_video_player.py` can be used.
## Testing
### Time window extraction
A small ground truth of speedups, slowdowns and lane change has been made. You can run unit tests by using the Python unittest module, for example:

Note: Each unit test uses a ground truth so it needs the respective original trip data to be in the 'preprocessing_input_directory' defined in the paths configuration.
```bash
python -m unittest -v acceleration_learning.tests.test_speedup_cutting
```
```bash
python -m unittest -v acceleration_learning.tests.test_slowdown_cutting
```
```bash
python -m unittest -v acceleration_learning.tests.test_window_maker
```
