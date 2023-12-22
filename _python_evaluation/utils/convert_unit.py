COEFFICIENT_TO_CONVERT_M_S_TO_KM_H = 3.6
KM_PER_TO_MILE_PER_H_COEF = 0.621371

def convert_km_per_h_to_m_per_s(speed_in_km_per_h: float) -> float:
    return speed_in_km_per_h / COEFFICIENT_TO_CONVERT_M_S_TO_KM_H

def convert_m_per_s_to_km_per_h(speed_in_m_per_s: float) -> float:
    return speed_in_m_per_s * COEFFICIENT_TO_CONVERT_M_S_TO_KM_H

def convert_km_per_h_to_mph(speed_in_km_per_h: float) -> float:
    return speed_in_km_per_h / KM_PER_TO_MILE_PER_H_COEF

def convert_m_per_s_to_mph(speed_in_m_per_s: float) -> float:
    return speed_in_m_per_s * 2.23694
