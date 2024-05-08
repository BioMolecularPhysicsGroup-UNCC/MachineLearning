#!/usr/bin/env python
# coding: utf-8

#%% 

#import tensorflow as tf
from tensorflow import keras
import numpy as np
import gc

#%%

################################# PARAMETERS #################################

training_species="Drosophila_2L_minimal"
Kmer_species=training_species
unknown_sequence = "A0A024G196_9STRA"
set_title = "DFT"
Neurons_per_layer = "10D"
Training_fraction = "p50"
window_size = 69
model_validation_number = 0 # enter which model to use if multiple models were trained for cross validation

########### CDS prediction (post processing) #################################

CDS_threshold_list = [3]
min_cds_length = 30
max_cds_length = 40

############################### END PARAMETERS ###############################

#%% 

###################### Data min-max scaler
# coded myself to work with arrays of the shape I am using

def Data_min_max_scaler(data_array):
    N = data_array.shape[0]
    nsensors = data_array.shape[1]
    scaled_array = np.zeros((N,nsensors))

    for j in range(0,nsensors):
        min0 = min(data_array[:,j])
        max0 = max(data_array[:,j])
        dif = max0 - min0
        if dif != 0:
            for k in range(0,N):
                x = ( data_array[k,j] - min0 ) / dif
                scaled_array[k,j] = x
        if dif == 0:
            for k in range(0,N):
                x = data_array[k,j]
                scaled_array[k,j] = x
    return scaled_array

#%%

def Adjust_predictions(minseq,maxseq,pred,seq,val_num):
    
    with open(seq,'r') as f:
        sequence = f.read()
        
    ### scan exon sites
    
    seq_len = len(sequence)
    
    potential_exon_sites = []
    
    for i in range(0,seq_len-max_cds_length):
        if sequence[i:i+3].upper() == "ATG":
            for j in range(i+min_cds_length,i+max_cds_length):
                if sequence[j:j+3].upper() in ["TAG","TAA","TGA"]:
                    if len(sequence[i:j+3]) % 3 == 0:
                        cds = sum(pred[i:j])
                        noncds = j - i - cds + 0.0001 - 1
                        potential_exon_sites.append((i,j,cds/noncds,"Single"))
    
                if sequence[j:j+2].upper() == "GT":
                    cds = sum(pred[i:j])
                    noncds = j - i - cds + 0.0001 - 1
                    potential_exon_sites.append((i,j,cds/noncds,"Start"))                   
    
        if sequence[i:i+2].upper() == "GT":
            for j in range(i+min_cds_length,i+max_cds_length):
                if sequence[j:j+3].upper() in ["TAG","TAA","TGA"]:
                    cds = sum(pred[i:j])
                    noncds = j - i - cds + 0.0001 - 1
                    potential_exon_sites.append((i,j,cds/noncds,"Stop"))
    
        if sequence[i:i+2].upper() == "AG":
            for j in range(i+min_cds_length,i+max_cds_length):
                if sequence[j:j+2].upper() == "GT":
                    cds = sum(pred[i:j])
                    noncds = j - i - cds + 0.0001 - 1
                    potential_exon_sites.append((i,j,cds/noncds,"Internal"))
                        
                        
    ### count potential exon sites backwards from largest to smallest
    ### Skip over small potential exons if larger predicted exon encompases it
    
    for threshold in CDS_threshold_list:
        
        predictions = np.zeros((seq_len))
        
        CDS_list = []
        
        exon_lengths = []
    
        single_exons = 0
        start_exons = 0
        stop_exons = 0
        internal_exons = 0
        
        total_cds = 0
        
        print(f"Minimum - Maximum cds length: {minseq}nt - {maxseq}nt")

        skip_start_site = 0
                        
        for k in range(len(potential_exon_sites)-1,0,-1):
            
            start = potential_exon_sites[k][0]
            stop = potential_exon_sites[k][1]
#            print(start,stop)
            
            if start == skip_start_site:
                continue
            
            if potential_exon_sites[k][2] > threshold:
                CDS_list.append((start,stop))
                skip_start_site = start
                for m in range(start,stop):
                    predictions[m] = 1
                if potential_exon_sites[k][3] == "Single":
                    single_exons += 1
                if potential_exon_sites[k][3] == "Start":
                    start_exons += 1
                if potential_exon_sites[k][3] == "Stop":
                    stop_exons += 1
                if potential_exon_sites[k][3] == "Internal":
                    internal_exons += 1
    
                  
        exon_start = 0                    
                   
        for i in range(0,seq_len-1):
            if predictions[i] == 0:
                if predictions[i+1] == 1:
                    exon_start = i+1
            if predictions[i] == 1:
                total_cds += 1
                if predictions[i+1] == 0:
                    exon_stop = i+1
                    exon_lengths.append(exon_stop-exon_start)
                
                
        cds_fraction = total_cds/seq_len
        
        mean_exon_length = np.mean(exon_lengths[:])
        
        print("Predicted CDS/nonCDS fraction at threshold =",threshold,":",cds_fraction)
        print("Mean exon length:",np.round(mean_exon_length,3))
        print("Single exons:",np.round(single_exons,3))
        print("Start exons:",np.round(start_exons,3))
        print("Stop exons:",np.round(stop_exons,3))  
        print("Internal exons:",np.round(internal_exons,3))
        print("--------------------------------------------------")
               
        hws = int( (window_size - 1 ) / 2 )
        
        with open(f"{unknown_sequence}_CDS_predictions.gff3", 'w') as f:

            for i in range(len(CDS_list)-1,0,-1):

                strand = "+"
                cds_start = CDS_list[i][0]
                cds_stop =  CDS_list[i][1]
                f.write(f"{set_title}\t")
                f.write("NN\t")
                f.write("CDS\t")
                f.write(f"{cds_start+hws}\t")
                f.write(f"{cds_stop+hws}\t.\t{strand}\t0")
                f.write(f'\tID=CDS:{i};Parent=transcript:FBtr0078171;protein_id=FBpp0077829\n')
                
                print(unknown_sequence,"CDS",cds_start+hws,cds_stop+hws)
                
        print("Output saved to:",f"{unknown_sequence}_CDS_predictions.gff3")
    
#%%

#results = []

def NN_predictions(name,val_num,model,x_tests,seq):
               
    pred_raw = model.predict(x_tests)
    pred = np.argmax(pred_raw, axis = 1)
    
    seq_len = pred.shape[0]    
    cds_count = 0
    for i in range(0,seq_len):
        if pred[i] == 1:
            cds_count += 1
    cds_fraction = cds_count/seq_len
    
    print("Initial predicted CDS/nonCDS fraction:",cds_fraction)
    
    ############ run CDS prediction algorithm
    
    Adjust_predictions(min_cds_length,max_cds_length,pred,seq,val_num)
    

#%%

model = keras.models.load_model(f"./Trained_models/Dense_model_{set_title}_{Neurons_per_layer}_{training_species}_{Training_fraction}_val_{model_validation_number}")

print("--------------------------------------------------")

#print("TEST PARAMETERS",training_species,Kmer_species,"ST",set_title,"NPL",Neurons_per_layer,"TF",Training_fraction,"WS",window_size)

test_data = np.load(f"./Test_datasets/{unknown_sequence}/{Kmer_species}/x_{unknown_sequence}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")

x_test = test_data[:,0:-1] # Sensor data
x_tests = Data_min_max_scaler(x_test)
del x_test
gc.collect()
    
seq = f"./Test_datasets/{unknown_sequence}/{Kmer_species}/s_{unknown_sequence}_{set_title}_ws{window_size}_ks_{Kmer_species}_sequence.txt"

NN_predictions(f"Dense_model_{set_title}_{Neurons_per_layer}_{training_species}_{Training_fraction}",
               model_validation_number,model,x_tests,seq)

print(f"END RESULTS {training_species}/{unknown_sequence} {Training_fraction}")

