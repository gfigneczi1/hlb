import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import pathlib
from gooey import GooeyParser, Gooey
from utils.utilities import utilities
from dynamics_utils.data_getter import data_getter
from dynamics_utils.dynamics_segmentation import dynamics_segmentation
from dynamics_utils.weighted_formulas import weighted_formulas

def save_weights(wf, formula):
    homedir = str(pathlib.Path(__file__).resolve().parent.parent)
    if not "segmented" in formula:
        out_list = [['w1', 'w2', 'w3', 'w4'],[wf.w1, wf.w2,wf.w3,wf.w4]]
    else:
        out_list = [['w11','w12','w13','w14','w21','w22','w23','w31','w32','w33','wInit', 'wMid', 'wSett'],
                        [wf.w11,wf.w12,wf.w13,wf.w14,wf.w21,wf.w22,wf.w23,wf.w31,wf.w32,wf.w33,wf.wInit, wf.wMid, wf.wSett]]
    df = pd.DataFrame(out_list)
    df.to_csv(os.path.join(homedir, "_temp", "Level_of_Dynamics","weights.csv"),index=False)
def process(path, wf, formula):
    #data getter is initialized
    dg = data_getter()
    a_y_rel, jerk = dg.get_continuous_data(path)
    #the lane change is segmented
    ds = dynamics_segmentation()
    df = pd.read_csv(path)
    t1, t2 = ds.segment_alternate(a_y_rel, df)
    t3 = len(a_y_rel)
    #empty initializations
    max_a_y_init, max_a_y_middle, max_a_y_sett = 0, 0, 0
    init_jerk, jerk_int_init, jerk_int_middle, jerk_int_sett = 0, 0, 0, 0
    init_time, middle_time, settling_time = 0, 0, 0
    #get initial phase data
    if t1 > 0:
        max_a_y_init = dg.get_max_a_y(a_y_rel, 0, t1)
        init_jerk = dg.get_init_jerk(jerk, 0, t1)
        jerk_int_init = dg.get_integral_jerk(jerk, 0, t1)
        init_time = t1
    #get middle phase data
    if t2 - t1 > 0:
        max_a_y_middle = dg.get_max_a_y(a_y_rel, t1, t2)
        jerk_int_middle = dg.get_integral_jerk(jerk, t1, t2)
        middle_time = t2 - t1
    #get settling phase data
    if t3 - t2 > 0:
        max_a_y_sett = dg.get_max_a_y(a_y_rel, t2, t3)
        jerk_int_sett = dg.get_integral_jerk(jerk, t2, t3)
    #get data spanning entire lane change
    max_a_y = max(max_a_y_init, max_a_y_middle, max_a_y_sett)
    jerk_int = jerk_int_init + jerk_int_middle + jerk_int_sett
    t11, t12, t13, t14 = ds.segment_ay_thr(a_y_rel) 
    t21, t22, t23, t24 = ds.segment_c01_thr(df)
    Tau = (t24) - (t12) #should be changed
    #the weighted formula that is used to calculate the Level of dynamics is different depending on the 'formula' input
    #all of the implementations of the weighted formula can be found in the weighted_formulas.py file
    if formula == 'weighted':
        LoD = wf.weighted_formula(max_a_y, init_jerk, jerk_int, Tau)
    elif formula == 'z_normalize':
        LoD = wf.z_normalize(max_a_y, init_jerk, jerk_int, Tau)
    elif formula == 'min_max_normalize':
        LoD = wf.min_max_normalize(max_a_y,  init_jerk, jerk_int, Tau)
    elif formula == 'segmented_weighted':
        LoD = wf.segmented_weighted(max_a_y_init, init_jerk, init_time, jerk_int_init,max_a_y_middle, jerk_int_middle, middle_time, max_a_y_sett, jerk_int_sett,settling_time)
    elif formula == 'segmented_z_normalize':
        LoD = wf.segmented_z_normalize(max_a_y_init, init_jerk, init_time, jerk_int_init,max_a_y_middle, jerk_int_middle, middle_time, max_a_y_sett, jerk_int_sett,settling_time)
    elif formula == 'segmented_min_max_normalize':
        LoD = wf.segmented_min_max_normalize(max_a_y_init, init_jerk, init_time, jerk_int_init,max_a_y_middle, jerk_int_middle,middle_time, max_a_y_sett, jerk_int_sett,settling_time)
    else:
        LoD = 0
    return LoD 
'''
This loops over all the files of a particular type of recording (natural, dynamic, comfort), calculates Level of Dynamics for all
of the files and creates a histogram with all of the results.
'''
def looping(wf, folder_path, RecordingType, formula):
    homedir = str(pathlib.Path(__file__).resolve().parent.parent)
    print(RecordingType)
    print(formula)
    all_LoDs = []
    for folder in os.listdir(folder_path):
        df = pd.DataFrame(["Level of Dynamics"])
        folder_LoDs = []
        for file in os.listdir(os.path.join(folder_path, folder)):
            path = os.path.join(folder_path, folder, file)
            LoD = process(path, wf, formula)
            all_LoDs.append(LoD)
            folder_LoDs.append(LoD)
            df.loc[len(df.index)] = [LoD]
        df.to_csv(os.path.join(homedir, "_temp", "Level_of_Dynamics",folder+'_LoD.csv'), index=False, header=False)
        plt.plot(folder_LoDs,color='blue')
        plt.title(folder + " Level of Dynamics development")
        plt.xlabel("Lane Changes")
        plt.ylabel("Level of Dynamics")
        plt.savefig(os.path.join(homedir, "_temp", "Level_of_Dynamics",folder+'_LoD_development.png'))
        plt.close()
    plt.hist(all_LoDs)
    print(np.average(all_LoDs))
    plt.title(RecordingType+" Level of Dynamics Distribution")
    plt.ylabel("# of Lane Changes")
    plt.xlabel("Level of Dynamics value")
    plt.savefig(os.path.join(homedir, "_temp", "Level_of_Dynamics",RecordingType+".png"))
    plt.show()

@Gooey(program_name='Weighted formula adjuster', default_size=(600, 650))
def main():
    parser = GooeyParser(description='input the values for the weighted formula')
    main_group = parser.add_argument_group("Main Options","Select main options")
    main_group.add_argument('--Formula', default='weighted', choices=['weighted','z_normalize','min_max_normalize','segmented_weighted', 'segmented_z_normalize', 'segmented_min_max_normalize'], widget='Dropdown')
    main_group.add_argument('--InputFolder', widget='DirChooser')
    simple_group = parser.add_argument_group("Weights for the simple formula","Adjust weights")
    simple_group.add_argument('--w1', default=1, widget='DecimalField')
    simple_group.add_argument('--w2', default=1, widget='DecimalField')
    simple_group.add_argument('--w3', default=1, widget='DecimalField')
    simple_group.add_argument('--w4', default=1, widget='DecimalField')
    segmented_group = parser.add_argument_group("Weights for segmented formula","Adjust weights")
    segmented_group.add_argument('--w11', default=1, widget='DecimalField')
    segmented_group.add_argument('--w12', default=1, widget='DecimalField')
    segmented_group.add_argument('--w13', default=1, widget='DecimalField')
    segmented_group.add_argument('--w14', default=1, widget='DecimalField')
    segmented_group.add_argument('--w21', default=1, widget='DecimalField')
    segmented_group.add_argument('--w22', default=1, widget='DecimalField')
    segmented_group.add_argument('--w23', default=1, widget='DecimalField')
    segmented_group.add_argument('--w31', default=1, widget='DecimalField')
    segmented_group.add_argument('--w32', default=1, widget='DecimalField')
    segmented_group.add_argument('--w33', default=1, widget='DecimalField')
    segmented_group.add_argument('--wInit', default=1, widget='DecimalField')
    segmented_group.add_argument('--wMid', default=1, widget='DecimalField')
    segmented_group.add_argument('--wSett', default=1, widget='DecimalField')
    args = parser.parse_args()
    folder_path = str(args.InputFolder)
    last_slash = 0
    for i, letter in enumerate(folder_path):
        if letter == '\\':
            last_slash = i
    RecordingType = folder_path[last_slash+1:]
    dg = data_getter()
    util = utilities()
    w1 = float(args.w1)
    w2 = float(args.w2)
    w3 = float(args.w3)
    w4 = float(args.w4)
    w11 = float(args.w11)
    w12 = float(args.w12)
    w13 = float(args.w13)
    w14 = float(args.w14)
    w21 = float(args.w21)
    w22 = float(args.w22)
    w23 = float(args.w23)
    w31 = float(args.w31)
    w32 = float(args.w32)
    w33 = float(args.w33)
    wInit = float(args.wInit)
    wMid = float(args.wMid)
    wSett = float(args.wSett)
    wf = weighted_formulas(w1, w2, w3, w4, w11, w12, w13, w14, w21, w22, w23,w31, w32, w33, wInit, wMid, wSett)
    if 'segmented' in args.Formula:
        print('w11=',w11,'w12=', w12,'w13=', w13, 'w14=', w14)
        print('w21=',w21,'w22=', w22, 'w23=', w23,)
        print('w31=',w31,'w32=', w32,'w33=', w33)
        print('w Init=',wInit,'w Mid=',wMid, 'w Sett=', wSett)
    else:
        print('w1=',w1,'w2=', w2,'w3=', w3, 'w4=', w4)
    homedir = str(pathlib.Path(__file__).resolve().parent.parent)
    path = os.path.join(homedir, "_temp", "Level_of_Dynamics")
    if not os.path.exists(path):
        os.makedirs(path)
    looping(wf, folder_path, RecordingType, args.Formula)
    save_weights(wf, args.Formula)


if __name__=="__main__":
    main()
            