#!/usr/bin/env python
# coding: utf-8

#%% 

import numpy as np
import os

#%%

################################# PARAMETERS #################################

species_list = ["Drosophila_2L_minimal"]
kmer_species_list = ["Drosophila_2L_minimal"]
window_size = 69
set_title = "DFT"

############################### END PARAMETERS ###############################

#%%

def Generate_test_data(x,test_sample_number,seq):
    
    sensor_data = x
    
    nsensors = sensor_data.shape[1]
    
    x_test = np.zeros((test_sample_number,nsensors))
    test_seq_list = []
    
    for i in range(0,test_sample_number):
        x_test[i] = sensor_data[i]
        test_seq_list.append(seq[i])
    
    test_seq = ''.join(test_seq_list)
            
    return x_test,test_seq

        
#%%
    
def Run_generate_test_dataset(species,ws,kmer_species,set_title):
    
    input_path = f"./Sensor_data/{species}/{kmer_species}/"
    
    parent_dir = "./"
    new_dir = "Test_datasets"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    parent_dir = "./Test_datasets"
    new_dir = f"{species}"    
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    parent_dir = f"./Test_datasets/{species}"
    new_dir = f"{kmer_species}"    
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    output_test_path = f"./Test_datasets/{species}/{kmer_species}/"
    
    sda = f"{input_path}{species}_ws{ws}_kmers_{kmer_species}_sensor_data.npy"
    seq = f"{input_path}{species}_ws{ws}_kmers_{kmer_species}_sensor_sequence.txt"
    
    sensor_data = np.load(sda)
    
    with open(seq,'r') as f:
        sensor_sequence = f.read()
    
    Sensor_data_samples = sensor_data.shape[0]
    
#%%

    x_test,test_seq = Generate_test_data(sensor_data,
                                                Sensor_data_samples,
                                                sensor_sequence)
    
#    print(x_test[0])
    
    output_test_filename = f"{species}_{set_title}_ws{ws}_ks_{kmer_species}"
    np.save(f"{output_test_path}x_{output_test_filename}_test",x_test)
    
    with open(f"{output_test_path}s_{output_test_filename}_sequence.txt", 'w') as f:
        f.write(test_seq)

#%%
        
for i in range(len(species_list)):
    for j in range(len(kmer_species_list)):        
        
        Run_generate_test_dataset(species_list[i],window_size,
                             kmer_species_list[j],set_title)
