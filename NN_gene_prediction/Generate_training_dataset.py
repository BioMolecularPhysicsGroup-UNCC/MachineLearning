#!/usr/bin/env python
# coding: utf-8

#%% 

import numpy as np
import random
import os

#%%

################################# PARAMETERS #################################

species_list = ["Drosophila_2L_minimal"]
window_size = 69
set_title = "DFT"
n_samples = 10
positive_fraction = 0.5

############################### END PARAMETERS ###############################
#%%
    
def Generate_training_data(x,y,tx,ty,start_train,stop_train):
    
    sensor_data = x
    sensor_labels = y
    
    targets = tx
    
    nsensors = sensor_data.shape[1]
    
    if sensor_data.shape[0] != sensor_labels.shape[0]:
        print("ERROR: Data and label arrays are not equal in length.")
    
    training_sample_number = stop_train - start_train
    
    x_train = np.zeros((training_sample_number,nsensors))
    y_train = np.zeros(training_sample_number)
    
    print(x_train.shape)
    
    for i in range(0,training_sample_number):
        target = int(targets[i])
        x_train[i] = sensor_data[target]
        y_train[i] = sensor_labels[target]
        
    y_train_column = y_train.reshape(-1,1)
        
    training_data_array = np.hstack((x_train, y_train_column))
            
    return training_data_array

       
#%%
    
def Run_generate_training_dataset(species,ws,set_title,
                         n_samples,positive_fraction):
    
    input_path = f"./Sensor_data/{species}/{species}/"
    
    Training_fraction = f"p{int(10*positive_fraction)}0"
    
    parent_dir = "./"
    new_dir = "Training_datasets"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    parent_dir = "./Training_datasets"
    new_dir = f"{species}"    
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    parent_dir = f"./Training_datasets/{species}"
    new_dir = f"{Training_fraction}"    
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    output_train_path = f"./Training_datasets/{species}/{Training_fraction}/"
    
    training_data_path = f"{input_path}{species}_ws{ws}_kmers_{species}_sensor_data.npy"
#    sla = f"{input_path}{chromosome}_ws{ws}_kmers_{kmer_species}_sensor_labels.npy"
    
    training_data_array = np.load(training_data_path)
    
    sensor_data = training_data_array[:,0:-1]
    sensor_labels = training_data_array[:,-1]
    
    Sensor_data_samples = sensor_data.shape[0]
    
    if sensor_data.shape[0] != sensor_labels.shape[0]:
        print("Data and label arrays are not equal in length.")
    
    ####################### Generate feature lists
    
    feature_list = []
    non_feature_list = []
    
    for i in range(0,Sensor_data_samples):
        if sensor_labels[i] == 0:
            non_feature_list.append(i)
        if sensor_labels[i] == 1:
            feature_list.append(i)
            
    ################# Geberate targets and labels
    
    Dataset_targets = np.zeros(n_samples)
    Dataset_labels = np.zeros(n_samples)
    
    if len(feature_list) > 0:
            
        for i in range(0,n_samples):
            pn = random.random()
            if pn < positive_fraction:
                target = random.choice(feature_list)
                Dataset_targets[i] = target
                Dataset_labels[i] = 1
                
            else:
                target = random.choice(non_feature_list)
                Dataset_targets[i] = target
                Dataset_labels[i] = 0
                
    if len(feature_list) == 0:
        print(f"EMPTY dataset {species}")

#%%
            
    training_set = Generate_training_data(sensor_data,sensor_labels,Dataset_targets,
                                             Dataset_labels,0,n_samples)
 
    output_train_filename = f"{species}_{Training_fraction}_ws{ws}_{set_title}"    
    np.save(f"{output_train_path}{output_train_filename}_training_set",training_set)

#%%

for i in range(len(species_list)):  
        
    Run_generate_training_dataset(species_list[i],window_size,set_title,
                                  n_samples,positive_fraction)
