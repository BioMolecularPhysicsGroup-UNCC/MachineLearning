#!/usr/bin/env python
# coding: utf-8

#%% 

import tensorflow as tf
from tensorflow import keras
import numpy as np
import gc
import os

#%%

################################# PARAMETERS #################################

Training_species="Drosophila_2L_minimal"
Kmer_species=Training_species
species_list = ["Drosophila_2L_minimal"]
set_title = "DFT"
Neurons_per_layer = "10D"
Training_fraction = "p50"
window_size = 69
cross_validation = 5

########### CDS prediction (post processing) #################################

initial_voting_threshold = 0.5
min_cds_length = 30
max_cds_length = 40
max_cds_threshold = 0.7
threshold_step = 0.2

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

def Adjust_predictions(minseq,maxseq,result_array,seq):
    
    with open(seq,'r') as f:
        sequence = f.read()
    
    y_test = result_array[0,:]
    pred = result_array[1,:]
    
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
                        if cds/noncds > initial_voting_threshold:
                            potential_exon_sites.append((i,j,cds/noncds,"Single"))
    
                if sequence[j:j+2].upper() == "GT":
                    cds = sum(pred[i:j])
                    noncds = j - i - cds + 0.0001 - 1
                    if cds/noncds > initial_voting_threshold:
                        potential_exon_sites.append((i,j,cds/noncds,"Start"))                   
    
        if sequence[i:i+3].upper() in ["AAG","CAG","TAG"]:
            for j in range(i+min_cds_length,i+max_cds_length):
                if sequence[j:j+3].upper() in ["TAG","TAA","TGA"]:
                    cds = sum(pred[i:j])
                    noncds = j - i - cds + 0.0001 - 1
                    if cds/noncds > initial_voting_threshold:
                        potential_exon_sites.append((i,j,cds/noncds,"Stop"))
    
        if sequence[i:i+3].upper() in ["AAG","CAG","TAG"]:
            for j in range(i+min_cds_length,i+max_cds_length):
                if sequence[j:j+2].upper() == "GT":
                    cds = sum(pred[i:j])
                    noncds = j - i - cds + 0.0001 - 1
                    if cds/noncds > initial_voting_threshold:
                        potential_exon_sites.append((i,j,cds/noncds,"Internal"))
                        
                        
    ### count potential exon sites backwards from largest to smallest
    ### Skip over small potential exons if larger predicted exon encompases it
    
    for threshold in np.arange(initial_voting_threshold,max_cds_threshold,threshold_step):
        
        predictions = np.zeros((seq_len))
        
        exon_lengths = []
    
        single_exons = 0
        start_exons = 0
        stop_exons = 0
        internal_exons = 0
        
        total_cds = 0
        
        print(f"Minimum - Maximum cds length: {minseq}nt - {maxseq}nt")
        print(f"CDS/nonCDS threshold: {threshold}")

        skip_start_site = 0
                        
        for k in range(len(potential_exon_sites)-1,0,-1):
            
            start = potential_exon_sites[k][0]
            stop = potential_exon_sites[k][1]
            
            if start == skip_start_site:
                continue
            
            if potential_exon_sites[k][2] > threshold:
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
        
        m = tf.keras.metrics.TruePositives()
        m.update_state(y_test, predictions)
        TP = m.result().numpy()
        
        m = tf.keras.metrics.FalsePositives()
        m.update_state(y_test, predictions)
        FP = m.result().numpy()
        
        m = tf.keras.metrics.TrueNegatives()
        m.update_state(y_test, predictions)
        TN = m.result().numpy()
        
        m = tf.keras.metrics.FalseNegatives()
        m.update_state(y_test, predictions)
        FN = m.result().numpy()
        
        
        # Sensitivity/recall = TP / ( TP + FN)
        
        sensitivity = TP / (TP + FN)
        print("Sensitivity: ",np.round(sensitivity,3))
        
        # Specicifity = TN / (TN + FP)
        
        specicifity = TN / (TN + FP)
        print("Specificity: ",np.round(specicifity,3))
        
        # F1 score = TP / ( TP + 0.5*(FP + FN) )
        
        F1 = TP / ( TP + FP/2 + FN/2 )
        print("F1 score: ",np.round(F1,3))
        
        print("Mean exon length:",np.round(mean_exon_length,3))
        print("Single exons:",np.round(single_exons,3))
        print("Start exons:",np.round(start_exons,3))
        print("Stop exons:",np.round(stop_exons,3))  
        print("Internal exons:",np.round(internal_exons,3))
        print("OUTPUT EXON",Training_species,Kmer_species,Test_species,"SN",np.round(sensitivity,3),"SP",np.round(specicifity,3),"F1",np.round(F1,3),"CDSF",np.round(cds_fraction,3),"THR",np.round(threshold,3))
        print("--------------------------------------------------")
    
#%%

#results = []

def test_nn(name,val_set,model,x_tests,y_test,seq):
        
    ################## Forward strand test data
    
    Accuracy = model.evaluate(x_tests, y_test, verbose=0)[1]
    
    pred = model.predict(x_tests)
    pred = np.argmax(pred, axis = 1)
    
    seq_len = pred.shape[0]    
    cds_count = 0
    for i in range(0,seq_len):
        if pred[i] == 1:
            cds_count += 1
    cds_fraction = cds_count/seq_len    
    
    result_array = np.concatenate((y_test,pred),axis=0)
    result_array.shape = (2,y_test.shape[0])
    
    parent_dir = "./"
    new_dir = "Dense_results"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    parent_dir = "./Dense_results"
    new_dir = f"{Training_species}"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    parent_dir = f"./Dense_results/{Training_species}"
    new_dir = f"{Test_species}"    
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
        
    parent_dir = f"./Dense_results/{Training_species}/{Test_species}"
    new_dir = f"{Training_fraction}"    
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
   
    ################### calculate TP, FP, TN and FN numbers
    
    m = tf.keras.metrics.TruePositives()
    m.update_state(y_test, pred)
    TP = m.result().numpy()
    
    m = tf.keras.metrics.FalsePositives()
    m.update_state(y_test, pred)
    FP = m.result().numpy()
    
    m = tf.keras.metrics.TrueNegatives()
    m.update_state(y_test, pred)
    TN = m.result().numpy()
    
    m = tf.keras.metrics.FalseNegatives()
    m.update_state(y_test, pred)
    FN = m.result().numpy()
    
    sensitivity = TP / (TP + FN)
    
    specicifity = TN / (TN + FP)
    
    F1 = TP / ( TP + FP/2 + FN/2 )
    
    print("OUTPUT NUCLEOTIDE",Training_species,Kmer_species,Test_species,"SN",np.round(sensitivity,3),"SP",np.round(specicifity,3),"F1",np.round(F1,3),"ACC",np.round(Accuracy,3),"CDSF",np.round(cds_fraction,3))
    
    ############ run CDS prediction algorithm
    
    Adjust_predictions(min_cds_length,max_cds_length,result_array,seq)
    
##%%
#
#Test_species="Drosophila_2L"
#models = []
#        
#for i in range(0,cross_validation):
#    models.append(keras.models.load_model(f"./Trained_models/Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}_validation_{i}"))    
#
#print("--------------------------------------------------")
#
#print("TEST PARAMETERS",Training_species,Kmer_species,"ST",set_title,"NPL",Neurons_per_layer,"TF",Training_fraction,"WS",window_size)
#
#x_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/x_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#x_tests = Data_min_max_scaler(x_test)
#del x_test
#gc.collect()
#
#y_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/y_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#
#seq = f"./Test_datasets/{Test_species}/{Kmer_species}/s_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_sequence.txt"
#
#for i in range(0,cross_validation):
#    test_nn(f"Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}",
#            i,models[i],x_tests,y_test,seq)
#   
#    
#print(f"END RESULTS {Training_species}/{Test_species} {Training_fraction}")
#
##%%
#
#Test_species="Human_21"
#models = []
#
#for i in range(0,cross_validation):
#    models.append(keras.models.load_model(f"./Trained_models/Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}_validation_{i}"))
#
#print("--------------------------------------------------")
#
#print("TEST PARAMETERS",Training_species,Kmer_species,"ST",set_title,"NPL",Neurons_per_layer,"TF",Training_fraction,"WS",window_size)
#
#x_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/x_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#x_tests = Data_min_max_scaler(x_test)
#del x_test
#gc.collect()
#
#y_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/y_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#
#seq = f"./Test_datasets/{Test_species}/{Kmer_species}/s_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_sequence.txt"
#
#for i in range(0,cross_validation):
#    test_nn(f"Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}",
#            i,models[i],x_tests,y_test,seq)
#
#print(f"END RESULTS {Training_species}/{Test_species} {Training_fraction}")
#
##%%
#
#Test_species="Mouse_19"
#models = []
#
#for i in range(0,cross_validation):
#    models.append(keras.models.load_model(f"./Trained_models/Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}_validation_{i}"))
#
#print("--------------------------------------------------")
#
#print("TEST PARAMETERS",Training_species,Kmer_species,"ST",set_title,"NPL",Neurons_per_layer,"TF",Training_fraction,"WS",window_size)
#
#x_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/x_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#x_tests = Data_min_max_scaler(x_test)
#del x_test
#gc.collect()
#
#y_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/y_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#
#seq = f"./Test_datasets/{Test_species}/{Kmer_species}/s_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_sequence.txt"
#
#for i in range(0,cross_validation):
#    test_nn(f"Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}",
#            i,models[i],x_tests,y_test,seq)
#
#print(f"END RESULTS {Training_species}/{Test_species} {Training_fraction}")
#
##%%
#
#Test_species="Worm_I"
#models = []
#
#for i in range(0,cross_validation):
#    models.append(keras.models.load_model(f"./Trained_models/Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}_validation_{i}"))
#
#print("--------------------------------------------------")
#
#print("TEST PARAMETERS",Training_species,Kmer_species,"ST",set_title,"NPL",Neurons_per_layer,"TF",Training_fraction,"WS",window_size)
#
#x_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/x_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#x_tests = Data_min_max_scaler(x_test)
#del x_test
#gc.collect()
#
#y_test = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/y_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")
#
#seq = f"./Test_datasets/{Test_species}/{Kmer_species}/s_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_sequence.txt"
#
#for i in range(0,cross_validation):
#    test_nn(f"Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}",
#            i,models[i],x_tests,y_test,seq)
#
#print(f"END RESULTS {Training_species}/{Test_species} {Training_fraction}")

#%%

for Test_species in species_list:
    
    models = []

    for i in range(0,cross_validation):
        models.append(keras.models.load_model(f"./Trained_models/Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}_val_{i}"))
    
    print("--------------------------------------------------")
    
    print("TEST PARAMETERS",Training_species,Kmer_species,"ST",set_title,"NPL",Neurons_per_layer,"TF",Training_fraction,"WS",window_size)
    
    test_data = np.load(f"./Test_datasets/{Test_species}/{Kmer_species}/x_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_test.npy")

    x_test = test_data[:,0:-1] # Sensor data
    x_tests = Data_min_max_scaler(x_test)
    del x_test
    gc.collect()
    
    y_test = test_data[:,-1] # labels
    
    seq = f"./Test_datasets/{Test_species}/{Kmer_species}/s_{Test_species}_{set_title}_ws{window_size}_ks_{Kmer_species}_sequence.txt"
    
    for i in range(0,cross_validation):
        test_nn(f"Dense_model_{set_title}_{Neurons_per_layer}_{Training_species}_{Training_fraction}",
                i,models[i],x_tests,y_test,seq)
    
    print(f"END RESULTS {Training_species}/{Test_species} {Training_fraction}")

