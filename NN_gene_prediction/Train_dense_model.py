#!/usr/bin/env python
# coding: utf-8

#%% 

import tensorflow as tf
import numpy as np
import gc
import os

#%%

################################# PARAMETERS #################################

species_list = ["Drosophila_2L_minimal"]
set_title = "DFT"
nepochs = 2
Neurons_per_layer = "10D"
Training_fraction = "p50"
window_size = 69
cross_validation = 5

############################### END PARAMETERS ###############################


#%% 

###################### Data min-max scaler

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

################################ Neural net

results = []

def train_nn(species,name,val_set):
    
    data_set = np.load(f"./Training_datasets/{species}/{Training_fraction}/{species}_{Training_fraction}_ws{window_size}_{set_title}_training_set.npy")
    
    x_train = data_set[:,0:-1]
    x_trains = Data_min_max_scaler(x_train)
    del x_train
    gc.collect()
    
    y_train = data_set[:,-1]

    val_set_size = int(x_trains.shape[0]/cross_validation)
    nsensors = x_trains.shape[1]
    
    model = tf.keras.models.Sequential([
      tf.keras.layers.Dense(10, input_shape=(nsensors,), activation='selu'),
      tf.keras.layers.Dense(2)])

    #model.summary()
    
    loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)
    
    model.compile(optimizer='adam',loss=loss_fn,
                  metrics=[tf.keras.metrics.SparseTopKCategoricalAccuracy(k=1)])

    val_start = val_set * val_set_size
    val_stop = (val_set + 1) * val_set_size
    
    x_trains_set = np.concatenate((x_trains[:val_start],x_trains[val_stop:]),axis=0)
    y_train_set = np.concatenate((y_train[:val_start],y_train[val_stop:]),axis=0)
    x_vals = x_trains[val_start:val_stop]
    y_val = y_train[val_start:val_stop]
    
    print("Training samples:",x_trains_set.shape)
    print("Validation samples:",x_vals.shape)
    
    for i in range(0,nepochs):
        model.fit(x_trains_set[:], y_train_set[:], epochs=1,
                  validation_data=(x_vals[:], y_val[:]),verbose=1)
        
    parent_dir = "./"
    new_dir = "Trained_models"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
        
    model.save(f"./Trained_models/{name}_val_{val_set}")
    
    del x_trains_set
    del y_train_set
    del x_vals
    del y_val
    del model
    gc.collect()
        
  
#%%
    
for i in range(len(species_list)):        
    for j in range(0,cross_validation):
        
        train_nn(species_list[i],f"Dense_model_{set_title}_10D_{species_list[i]}_{Training_fraction}",j)