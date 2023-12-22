from dynamics_utils.norm_params_getter import norm_params_getter
class weighted_formulas():
    def __init__(self, w1, w2, w3, w4, w11, w12, w13, w14, w21, w22,w23, w31, w32, w33, wInit, wMid, wSett):
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.w4 = w4
        self.w11 = w11
        self.w12 = w12
        self.w13 = w13
        self.w14 = w14
        self.w21 = w21
        self.w22 = w22
        self.w23 = w23
        self.w31 = w31
        self.w32 = w32
        self.w33 = w33
        self.wInit = wInit
        self.wMid = wMid
        self.wSett = wSett
    def z_norm(self,value, mu, std):
        norm = (value-mu)/std
        return norm
    def min_max_norm(self,value, min, max):
        norm = (value-min)/(max-min)
        return norm
    def weighted_formula(self, max_a_y, init_jerk,jerk_int, Tau):
        LoD = self.w1 * init_jerk + self.w2 * max_a_y + self.w3 * jerk_int 
        if Tau > 0:
            LoD += self.w4*(1/Tau)
        return LoD
    def z_normalize(self, max_a_y, init_jerk, jerk_int, Tau):
        npg = norm_params_getter()
        a_mu, a_std = npg.get_z_norm_params('max_a_y')
        #overall int jerk
        Ji_mu, Ji_std = npg.get_z_norm_params('jerk_int')
        #lane change time
        T_mu, T_std = npg.get_z_norm_params('Tau')
        #init jerk
        mu_init_jerk, std_init_jerk = npg.get_z_norm_params('init_jerk')
        LoD = self.w1 * self.z_norm(init_jerk, mu_init_jerk, std_init_jerk) + self.w2 * self.z_norm(max_a_y, a_mu, a_std) + self.w3 * self.z_norm(jerk_int, Ji_mu, Ji_std) 
        if Tau > 0:
            LoD += self.w4*self.z_norm(1/Tau,T_mu, T_std)
        return LoD
    def min_max_normalize(self, max_a_y, init_jerk, jerk_int, Tau):
        npg = norm_params_getter()
        #overall max acceleration
        a_min, a_max = npg.get_mm_norm_params('max_a_y')
        #overall int jerk
        Ji_min, Ji_max = npg.get_mm_norm_params('jerk_int')
        #lane change time
        T_min, T_max = npg.get_mm_norm_params('Tau')
        #init jerk
        min_init_jerk, max_init_jerk = npg.get_mm_norm_params('init_jerk')
        LoD = self.w1 * self.min_max_norm(init_jerk, min_init_jerk, max_init_jerk) + self.w2 * self.min_max_norm(max_a_y, a_min, a_max) + self.w3 * self.min_max_norm(jerk_int, Ji_min, Ji_max) 
        if Tau > 0:
            LoD += self.w4*self.min_max_norm(1/Tau,T_min, T_max)
        return LoD
    def segmented_weighted(self, max_a_y_init, init_jerk, init_time,jerk_int_init, max_a_y_middle, jerk_int_middle,middle_time, max_a_y_sett, jerk_int_sett,settling_time):
        init_LoD = self.w11 * max_a_y_init
        init_LoD = init_LoD + self.w12 * init_jerk
        init_LoD = init_LoD + self.w13 * init_time
        init_LoD = init_LoD + self.w14 * jerk_int_init
        middle_LoD = self.w21 * max_a_y_middle
        middle_LoD = middle_LoD + self.w22 * jerk_int_middle
        middle_LoD = middle_LoD + self.w23 * middle_time
        sett_LoD = self.w31 * max_a_y_sett 
        sett_LoD= sett_LoD + self.w32 * jerk_int_sett 
        sett_LoD = sett_LoD + self.w33* settling_time
        LoD = (self.wInit*init_LoD + self.wMid*middle_LoD +self.wSett*sett_LoD)
        return LoD
    def segmented_z_normalize(self, max_a_y_init, init_jerk, init_time, jerk_int_init, max_a_y_middle, jerk_int_middle,middle_time, max_a_y_sett, jerk_int_sett,settling_time):
        npg = norm_params_getter()
        mu_max_a_y_init , std_max_a_y_init = npg.get_z_norm_params('max_a_y_init')
        #
        mu_init_jerk, std_init_jerk = npg.get_z_norm_params('init_jerk')
        #
        mu_init_time, std_init_time = npg.get_z_norm_params('init_tau')
        #
        mu_jerk_int_init, std_jerk_int_init = npg.get_z_norm_params('jerk_int_init')
        '''MIDDLE'''
        mu_max_a_y_middle , std_max_a_y_middle = npg.get_z_norm_params('max_a_y_middle')
        #
        mu_jerk_int_middle , std_jerk_int_middle = npg.get_z_norm_params('jerk_int_middle')
        #
        mu_mid_time, std_mid_time = npg.get_z_norm_params('middle_tau')
        '''SETTLING'''
        mu_max_a_y_sett , std_max_a_y_sett = npg.get_z_norm_params('max_a_y_sett')
        #
        mu_jerk_int_sett , std_jerk_int_sett = npg.get_z_norm_params('jerk_int_sett')
        ''''WE Should replace it with 1/T data'''
        mu_sett_time , std_sett_time = npg.get_z_norm_params('settling_tau')
        #
        init_LoD = self.w11 * self.z_norm(init_jerk, mu_init_jerk, std_init_jerk) 
        init_LoD = init_LoD + self.w12 * self.z_norm(max_a_y_init, mu_max_a_y_init, std_max_a_y_init) 
        init_LoD = init_LoD + self.w13 * self.z_norm(init_time, mu_init_time, std_init_time)
        init_LoD = init_LoD + self.w14 * self.z_norm(jerk_int_init, mu_jerk_int_init, std_jerk_int_init)
        middle_LoD = self.w21 * self.z_norm(max_a_y_middle, mu_max_a_y_middle, std_max_a_y_middle) 
        middle_LoD= middle_LoD + self.w22 * self.z_norm(jerk_int_middle, mu_jerk_int_middle, std_jerk_int_middle)
        middle_LoD = middle_LoD + self.w23 * self.z_norm(middle_time, mu_mid_time, std_mid_time)
        sett_LoD = self.w31 * self.z_norm(max_a_y_sett, mu_max_a_y_sett, std_max_a_y_sett) 
        sett_LoD= sett_LoD + self.w32 * self.z_norm(jerk_int_sett, mu_jerk_int_sett, std_jerk_int_sett) 
        sett_LoD = sett_LoD + self.w33* self.z_norm(settling_time, mu_sett_time, std_sett_time)
        LoD = (self.wInit*init_LoD + self.wMid*middle_LoD +self.wSett*sett_LoD)
        return LoD
    def segmented_min_max_normalize(self, max_a_y_init, init_jerk, init_time,jerk_int_init, max_a_y_middle, jerk_int_middle,middle_time, max_a_y_sett, jerk_int_sett,settling_time):
        npg = norm_params_getter()
        min_max_a_y_init, max_max_a_y_init = npg.get_mm_norm_params('max_a_y_init')
        #
        min_init_jerk, max_init_jerk = npg.get_mm_norm_params('init_jerk')
        #
        min_init_time, max_init_time = npg.get_mm_norm_params('init_tau')
        #
        min_jerk_int_init, max_jerk_int_init = npg.get_mm_norm_params('jerk_int_init')
        '''MIDDLE'''
        min_max_a_y_middle, max_max_a_y_middle = npg.get_mm_norm_params('max_a_y_middle')
        #
        min_jerk_int_middle, max_jerk_int_middle = npg.get_mm_norm_params('jerk_int_middle')
        #
        min_mid_time, max_mid_time = npg.get_mm_norm_params('middle_tau')
        '''SETTLING'''
        min_max_a_y_sett, max_max_a_y_sett = npg.get_mm_norm_params('max_a_y_sett')
        #
        min_jerk_int_sett, max_jerk_int_sett = npg.get_mm_norm_params('jerk_int_sett')
        #
        min_sett_time, max_sett_time = npg.get_mm_norm_params('settling_tau')
        #
        init_LoD = self.w11 * self.min_max_norm(max_a_y_init, min_max_a_y_init, max_max_a_y_init) 
        init_LoD = init_LoD + self.w12 * self.min_max_norm(init_jerk, min_init_jerk, max_init_jerk) 
        init_LoD = init_LoD + self.w13 * self.min_max_norm(init_time, min_init_time, max_init_time)
        init_LoD = init_LoD + self.w14 * self.min_max_norm(jerk_int_init, min_jerk_int_init, max_jerk_int_init)
        middle_LoD = self.w21 * self.min_max_norm(max_a_y_middle, min_max_a_y_middle, max_max_a_y_middle) 
        middle_LoD = middle_LoD + self.w22 * self.min_max_norm(jerk_int_middle, min_jerk_int_middle, max_jerk_int_middle)
        middle_LoD = middle_LoD + self.w23 * self.min_max_norm(middle_time, min_mid_time, max_mid_time)
        sett_LoD = self.w31 * self.min_max_norm(max_a_y_sett, min_max_a_y_sett, max_max_a_y_sett) 
        sett_LoD= sett_LoD + self.w32 * self.min_max_norm(jerk_int_sett, min_jerk_int_sett, max_jerk_int_sett) 
        sett_LoD = sett_LoD + self.w33* self.min_max_norm(settling_time, min_sett_time, max_sett_time)
        LoD = (self.wInit*init_LoD + self.wMid*middle_LoD +self.wSett*sett_LoD)
        return LoD