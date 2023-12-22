from enum import Enum

class DrivingEventEnum(Enum):
    UPSTART = "upstart"
    STOP = "stop"
    TURN = "turn"
    SPEED_UP = "speed_up"
    SLOW_DOWN = "slow_down"
    CRUISE = "cruise"
    LANE_CHANGE = "lane_change"
