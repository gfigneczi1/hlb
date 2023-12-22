from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
build_options = {'packages': [], 'excludes': [], 'includes': ["conversion/reader/read_mdf.py"]}

import sys
base = 'Win32GUI' if sys.platform=='win32' else None

executables = [
    Executable('main.py', base=base, target_name = 'conv4kpi')
]

setup(name='conv4kpi',
      version = '1.0',
      description = 'Converter for the online KPI tool',
      options = {'build_exe': build_options},
      executables = executables)
