from launch import LaunchDescription
from launch_ros.actions import Node

def generate_launch_description():
    return LaunchDescription([
        Node(
            package='do_mpc_ros2',
            executable='vehicle_simu',
            name='mpc_vehicle_simu_node',
            output='screen'
        ),
        Node(
            package='do_mpc_ros2',
            executable='mpc_controller',
            name='mpc_node',
            output='screen'
        )
    ])
