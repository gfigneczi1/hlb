from setuptools import find_packages, setup
import os
import glob

package_name = 'do_mpc_ros2'

setup(
    name=package_name,
    version='0.0.0',
    packages=find_packages(exclude=['test']),
    data_files=[
        ('share/ament_index/resource_index/packages',
            ['resource/' + package_name]),
        ('share/' + package_name, ['package.xml']),
        ('share/' + package_name, ['launch/launch_mpc_sim.py']),
    ],
    install_requires=['setuptools'],
    zip_safe=True,
    maintainer='dobayt',
    maintainer_email='dobayt@todo.todo',
    description='TODO: Package description',
    license='Apache-2.0',
    tests_require=['pytest'],
    entry_points={
        'console_scripts': [
            'mpc_controller = do_mpc_ros2.run_controller:main',
            'vehicle_simu = do_mpc_ros2.simu_model:main'
        ],
    },
)
