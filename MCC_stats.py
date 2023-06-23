import sys
import numpy as np
from re import finditer
import matplotlib.pyplot as plt
from tabulate import tabulate




def calc_MCC_F_val(prediction, actual):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    prediction_bp=set()
    reference_bp=set()
    
    for i, (a, b) in enumerate(zip(prediction, actual)):
            if a != -1:
                if i<a:
                    prediction_bp.add((i,a))
                else:
                    prediction_bp.add((a,i))
            if b != -1:
                if i<b:
                    reference_bp.add((i,b))
                else:
                    reference_bp.add((b,i))
    for (i,j) in reference_bp:
        if (i,j) in prediction_bp:
            TP += 1
        else:
            FN+=1
    for i in range(0,len(actual)-1):
        for j in range (i+1, len(actual)):
            if (i,j) not in reference_bp and (i,j) not in prediction_bp:
                TN += 1

    for (i,j) in prediction_bp:
        if (i,j) not in reference_bp:
            FP+=1

    PPV = TP / (TP + FP)
    Sensitivity = TP / (TP + FN)
    MCC = ((TP * TN) - (FP * FN)) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    F_Value = 0.5 * (PPV + Sensitivity)

    return MCC, F_Value, PPV, Sensitivity


def parse_dot_bracket(input):
    output = np.full(len(input), -1)
    more = True
    while more:
        more = False

        #find matched parenthesis
        for x in finditer(r"\([^()]*\)", input):
            more = True
            output[x.start()] = x.end()-1
            output[x.end()-1] = x.start()

            input=input[0:x.start()] + "." + input[x.start()+1:x.end()-1] + "." + input[x.end():]

    return output


def calc_overall_MCC_F_val_PPV_Sensitivity(predictions, references):
    PPVs=[]
    Sensitivities=[]
    F_vals=[]
    MCCs=[]   
    for transcript_name in references.keys():
        TP = 0
        TN = 0
        FP = 0
        FN = 0
        prediction_bp=set()
        reference_bp=set()
        reference = parse_dot_bracket(references[transcript_name][1])
        prediction = parse_dot_bracket(predictions[transcript_name][1]) 
    
    
        for i, (a, b) in enumerate(zip(prediction, reference)):
            if a != -1:
                if i<a:
                    prediction_bp.add((i,a))
                else:
                    prediction_bp.add((a,i))
            if b != -1:
                if i<b:
                    reference_bp.add((i,b))
                else:
                    reference_bp.add((b,i))
        
        for (i,j) in reference_bp:
            if (i,j) in prediction_bp:
                TP += 1
            else:
                FN+=1
        for i in range(0,len(reference)-1):
            for j in range (i+1, len(reference)):
                if (i,j) not in reference_bp and (i,j) not in prediction_bp:
                    TN += 1
        for (i,j) in prediction_bp:
            if (i,j) not in reference_bp:
                FP+=1

        PPV = TP / (TP + FP)
        Sensitivity = TP / (TP + FN)
        MCC = ((TP * TN) - (FP * FN)) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        F_Value = 0.5 * (PPV + Sensitivity)

        PPVs.append(PPV)
        Sensitivities.append(Sensitivity)
        MCCs.append(MCC)
        F_vals.append(F_Value)

    MCC_all=sum(MCCs)/len(references.keys())
    F_Value_all=sum(F_vals)/len(references.keys())
    PPV_all=sum(PPVs)/len(references.keys())
    Sensitivity_all=sum(Sensitivities)/len(references.keys())
    
    return MCC_all, F_Value_all, PPV_all, Sensitivity_all, MCCs

def scatterplot(prediction_1, prediction_2, condition1, condition2):
    #scatterplot MCC
  
    plt.rcParams.update({'font.size': 16})
    plt.rcParams.update({'axes.titlesize': 16})
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    fig = plt.figure(figsize=(6,6))
    ax = fig.gca()
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="#0571b0")
    for transcript_name in prediction_1.keys():
        ax.scatter(prediction_1[transcript_name][2], prediction_2[transcript_name][2], color="#ca0020", alpha=0.7)
    ax.set_xlim(left=0.0, right=1.0)
    ax.set_ylim(bottom=0.0, top=1.0)
    # ax.set_xlabel(condition1)
    # ax.set_ylabel(condition2)
    ax.set_xlabel("MCC unstrained")
    ax.set_ylabel("MCC constrained")
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    # plt.legend(by_label.values(), by_label.keys())
    #plt.legend()
    plt.tight_layout()
    # plt.show()
    plt.savefig("./MCC_scatter.pdf", dpi=600)
    plt.close()

def stats_table(condition1, MCC_prediction_1_all, F_val_prediction_1_all, PPV_prediction_1_all, Sens_prediction_1_all, 
    condition2, MCC_prediction_2_all, F_val_prediction_2_all, PPV_prediction_2_all, Sens_prediction_2_all):
    
    stats=[[condition1, MCC_prediction_1_all, F_val_prediction_1_all, PPV_prediction_1_all, Sens_prediction_1_all], 
    [condition2, MCC_prediction_2_all, F_val_prediction_2_all, PPV_prediction_2_all, Sens_prediction_2_all]]

    print(tabulate(stats, headers=["prediction", "MCC", "F-val", "PPV", "Sensitivity"]))
    print("\n")

def read_file(inputfile):
    dict = {}
    line_counter=0
    current_id = ""
    with open(inputfile) as file:
        for line in file:
            line_counter += 1
            line=line.rstrip()
            if line.startswith(">"):        
                name=line[1:]               #transcript_name
                line_counter = 1
                dict[name]=[]
                current_id=name
            if line_counter == 2:           # RNA sequence
                dict[current_id].append(line)
            if line_counter == 3:           # structure
                dict[current_id].append(line)
        line_counter = 0
        current_id = ""
    return dict


if __name__ == "__main__":
    
    prediction_condition1_file = sys.argv[1]
    condition1 = sys.argv[2]
    prediction_condition2_file = sys.argv[3]
    condition2 = sys.argv[4]
    ref_structures_file = sys.argv[5]

    """store the structures in dicts"""

    prediction_1 = read_file(prediction_condition1_file)
    prediction_2 = read_file(prediction_condition2_file)
    ref_structures = read_file(ref_structures_file)


    improved_folding=[]
    for transcript_name in prediction_1.keys():
        reference = parse_dot_bracket(ref_structures[transcript_name][1])
        prediction_1_parsed = parse_dot_bracket(prediction_1[transcript_name][1])
        prediction_2_parsed = parse_dot_bracket(prediction_2[transcript_name][1])
        MCC_prediction_1, F_prediction_1,_ , _ = calc_MCC_F_val(prediction_1_parsed, reference)
        MCC_prediction_2, F_prediction_2,_ , _ = calc_MCC_F_val(prediction_2_parsed, reference)

 
        if MCC_prediction_1 <= MCC_prediction_2:
            improved_folding.append(transcript_name)

        prediction_1[transcript_name].append(MCC_prediction_1)
        prediction_1[transcript_name].append(F_prediction_1)
        prediction_2[transcript_name].append(MCC_prediction_2)
        prediction_2[transcript_name].append(F_prediction_2)

    scatterplot(prediction_1, prediction_2, condition1, condition2)

    # print out overall MCC, F-val, PPV, Sensitivity
    MCC_prediction_1_all, F_val_prediction_1_all, PPV_prediction_1_all, Sens_prediction_1_all,_ = calc_overall_MCC_F_val_PPV_Sensitivity(prediction_1, ref_structures)
    MCC_prediction_2_all, F_val_prediction_2_all, PPV_prediction_2_all, Sens_prediction_2_all,_ = calc_overall_MCC_F_val_PPV_Sensitivity(prediction_2, ref_structures)

    stats_table(condition1, MCC_prediction_1_all, F_val_prediction_1_all, PPV_prediction_1_all, Sens_prediction_1_all, 
    condition2, MCC_prediction_2_all, F_val_prediction_2_all, PPV_prediction_2_all, Sens_prediction_2_all)

    print("\n")
    print("improved structure prediction : "+str(improved_folding))
