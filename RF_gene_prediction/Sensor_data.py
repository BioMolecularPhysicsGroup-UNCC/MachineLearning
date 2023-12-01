#!/usr/bin/env python
# coding: utf-8

#%% 

import numpy as np
from collections import Counter
#import random
import os
import pickle

#%%

################################# USER PARAMETERS #################################

### For species parameter, enter the name of the target fasta file wihthout file extension

species = "Drosophila_2L_minimal"

################################# OTHER PARAMETERS #################################

kmer_species_list = ["Drosophila_2L","Human_21","Worm_I",
                     "Arabidopsis_1","Saccharomyces_IV"]
window_size = 99
n_sensors = 14

############################## Sensor functions ###############################
# All sensors are functions that take strings of letters (A,C,G,T) as input 
# and output numbers which are the features used by the model.

# A content sensor
# seq: input sequence in string format
# outputs A content of input sequence

def A_content_sensor(seq):
    A_count = 0
    for i in seq:
        if i.upper() == 'A':
            A_count +=1
    return A_count/len(seq)

# C content sensor
# seq: input sequence in string format
# outputs C content of input sequence

def C_content_sensor(seq):
    C_count = 0
    for i in seq:
        if i.upper() == 'C':
            C_count +=1
    return C_count/len(seq)

# D content sensor
# seq: input sequence in string format
# outputs G content of input sequence

def G_content_sensor(seq):
    G_count = 0
    for i in seq:
        if i.upper() == 'A':
            G_count +=1
    return G_count/len(seq)

# T content sensor
# seq: input sequence in string format
# outputs T content of input sequence

def T_content_sensor(seq):
    T_count = 0
    for i in seq:
        if i.upper() == 'A':
            T_count +=1
    return T_count/len(seq)

################# A autocorrelation sensor
# seq: input sequence in string format
# outputs max/min autocorrelation value
    
def A_autocorrelation_sensor(seq):
    A1_count = 0
    A2_count = 0
    A3_count = 0
    L = len(seq)
    for i in range(0,L-1,3):
        if seq[i].upper() == 'A':
            A1_count +=1
        if seq[i+1].upper() == 'A':
            A2_count +=1
        if seq[i+2].upper() == 'A':
            A3_count +=1
            
    numerator = max(A1_count,A2_count,A3_count)
    denominator = min(A1_count,A2_count,A3_count) + 1
            
    return numerator/denominator

# C autocorrelation sensor
# seq: input sequence in string format
# outputs max/min autocorrelation value
    
def C_autocorrelation_sensor(seq):
    C1_count = 0
    C2_count = 0
    C3_count = 0
    L = len(seq)
    for i in range(0,L-1,3):
        if seq[i].upper() == 'C':
            C1_count +=1
        if seq[i+1].upper() == 'C':
            C2_count +=1
        if seq[i+2].upper() == 'C':
            C3_count +=1
            
    numerator = max(C1_count,C2_count,C3_count)
    denominator = min(C1_count,C2_count,C3_count) + 1
            
    return numerator/denominator

# G autocorrelation sensor
# seq: input sequence in string format
# outputs max/min autocorrelation value
    
def G_autocorrelation_sensor(seq):
    G1_count = 0
    G2_count = 0
    G3_count = 0
    L = len(seq)
    for i in range(0,L-1,3):
        if seq[i].upper() == 'G':
            G1_count +=1
        if seq[i+1].upper() == 'G':
            G2_count +=1
        if seq[i+2].upper() == 'G':
            G3_count +=1
            
    numerator = max(G1_count,G2_count,G3_count)
    denominator = min(G1_count,G2_count,G3_count) + 1
            
    return numerator/denominator

# T autocorrelation sensor
# seq: input sequence in string format
# outputs max/min autocorrelation value
    
def T_autocorrelation_sensor(seq):
    T1_count = 0
    T2_count = 0
    T3_count = 0
    L = len(seq)
    for i in range(0,L-1,3):
        if seq[i].upper() == 'T':
            T1_count +=1
        if seq[i+1].upper() == 'T':
            T2_count +=1
        if seq[i+2].upper() == 'T':
            T3_count +=1
            
    numerator = max(T1_count,T2_count,T3_count)
    denominator = min(T1_count,T2_count,T3_count) + 1
            
    return numerator/denominator


# Kmer sensor
    
class Kmer_sensor():
    
    def __init__(self, coding_freqs=None, noncoding_freqs=None,k=4):
        
        self.coding_freqs = Counter()
        self.noncoding_freqs = Counter()
        self.tuple_size = k

                
    def encode(self, seq):
        # Get all words
        seqstr = ""
        seqstr = seqstr.join(seq).upper()
        words = self._count_kmers(seqstr, self.tuple_size)

        # For each word
        pref_sum = 0
        for word in dict(words).keys():
            cfreq  = self.coding_freqs.get(word, 0)
            ncfreq = self.noncoding_freqs.get(word, 0)

            if cfreq == 0 or ncfreq == 0: continue
            pref_sum += np.log(cfreq / ncfreq)

        # Return sum of coding preference values
        return pref_sum

    def _count_kmers(self,seq, k0):
        """Count kmers in sequence"""
        return Counter(seq[start:start+k0] for start in range(len(seq) - k0 + 1))
    
# Shannon entropy sensor

def Entropy_sensor(seq):
    nA = 0
    nT = 0
    nC = 0
    nG = 0
    l = len(seq)
    for i in range(0,l):
        if seq[i].upper() == 'A':
            nA += 1
        if seq[i].upper() == 'T':
            nT += 1
        if seq[i].upper() == 'C':
            nC += 1
        if seq[i].upper() == 'G':
            nG += 1

    pA = nA/l
    pT = nT/l
    pC = nC/l
    pG = nG/l

    E = 0

    if pA > 0:
        E += -pA * np.log(pA)
    if pT > 0:
        E += -pT * np.log(pT)
    if pC > 0:
        E += -pC * np.log(pC)
    if pG > 0:
        E += -pG * np.log(pG)
    
    SE = E*l
    
    return SE

#%% 

############################## Sensor sweep code ##############################
# Outputs array of sensor output vectors
# cdi = dna start coordinate
# cdf = dna stop coordinate
# seq: input dna sequence from fasta
# ws: window size

def Sensor_sweep(seq, target_nt,ws,k6_sensor,k5_sensor,k4_sensor,k3_sensor,k2_sensor,nsensors):
    output_array = np.zeros(nsensors)
    r = int((ws-1)/2) # half window size
    row = 0
    i = target_nt
    target_sequence = seq[i-r:i+r+1]
    output_array[row] = A_content_sensor(target_sequence); row = row + 1;
    output_array[row] = C_content_sensor(target_sequence); row = row + 1;
    output_array[row] = G_content_sensor(target_sequence); row = row + 1;
    output_array[row] = T_content_sensor(target_sequence); row = row + 1;
    output_array[row] = k6_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k5_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k4_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k3_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k2_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = Entropy_sensor(target_sequence); row = row + 1;
    output_array[row] = A_autocorrelation_sensor(target_sequence); row = row + 1;
    output_array[row] = C_autocorrelation_sensor(target_sequence); row = row + 1;
    output_array[row] = G_autocorrelation_sensor(target_sequence); row = row + 1;
    output_array[row] = T_autocorrelation_sensor(target_sequence); row = row + 1;
    return output_array

#%%
    
############################## Data generator ##############################

# This class is used to prepare datasets for the neural network
# fasta file and gff file must be in the same directory as this code
# Default values of class variables are set to use a small test dataset
    
class Data_generator():
    
    def __init__(self,ws,ns):
        
        self.fasta_input = "Drosophila_2L_minimal.fna" # my filename for the input fasta file
#        self.gff_input = "Drosophila_2L_minimal.gff3" # my filename for the input gff file
#        self.kmer_fasta_input = "Drosophila_2L_minimal.fna" # my filename for the input fasta file
#        self.kmer_gff_input = "Drosophila_2L_minimal.gff3" # my filename for the input gff file
        
        self.feature = "CDS" # sets type of feature to consider for positive samples
        self.nsensors = ns
#        self.chromosome = "2L"
#        self.kmer_chromosome = "2L"
        self.window_size = ws
#        self.kmer_training_multiplier = ktm
        
# Reads a gff file and stores the start and stop coordinates for each feature of the given type
# Reads a fasta file to generate a DNA string object for use by other functions
# Generates a label string the same length as the DNA string object
# Label string has '1' where there is a feature of the given type (eg. CDS)
# Label string has '0' where there is no feature of the given type (eg. CDS)
# Locations of features (eg. CDS) are extracted from inout gff file
        
    def Generate_sequence_and_labels(self):

                               
        genome = open(self.fasta_input)

        lines = genome.readlines()
        DNASeq0 = []
        for i in range (0, len(lines)):
            if lines[i][0:1] != ">":
                DNASeq0.append(lines[i].strip("\n"))
        genome.close()
        
        self.DNASeq = ''.join(DNASeq0)
        
                    
############# Function to load kmer sensors

    def load_kmer_sensors(self,ks):
        
        with open('./Kmer_sensors/kmer_6_sensor_' + ks + '.pkl', 'rb') as f:
            self.kmer_6_sensor = pickle.load(f)

        with open('./Kmer_sensors/kmer_5_sensor_' + ks + '.pkl', 'rb') as f:
            self.kmer_5_sensor = pickle.load(f)
            
        with open('./Kmer_sensors/kmer_4_sensor_' + ks + '.pkl', 'rb') as f:
            self.kmer_4_sensor = pickle.load(f)
        
        with open('./Kmer_sensors/kmer_3_sensor_' + ks + '.pkl', 'rb') as f:
            self.kmer_3_sensor = pickle.load(f)
        
        with open('./Kmer_sensors/kmer_2_sensor_' + ks + '.pkl', 'rb') as f:
            self.kmer_2_sensor = pickle.load(f)

# Produces sensor output data over a specific range of nucleotides
# Also produces a labeled string for comparison to predictions
    
    def Run_sensors(self,startnt,stopnt):
        
       
        test_seq = self.DNASeq[startnt:stopnt]
        
        test_number = len(test_seq)
                
        x_test0 = np.zeros((test_number,self.nsensors)) # all sensors

        for i in range(0,test_number):
            target = startnt + i
            x_test0[i] = Sensor_sweep(self.DNASeq, target,self.window_size,
                                      self.kmer_6_sensor,self.kmer_5_sensor,self.kmer_4_sensor,
                                      self.kmer_3_sensor,self.kmer_2_sensor,self.nsensors)
                
        return x_test0,test_seq      
    
#%%

def Generate_sensor_data(species,window_size,nsensors,kmer_species):
    
    hws = int( (window_size - 1 ) / 2 )

    Data = Data_generator(window_size,nsensors)
    
    parent_dir = "./"
    new_dir = "Sensor_data"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
        
    Data.fasta_input = species + ".fna"
    
    Data.load_kmer_sensors(kmer_species)
  
    parent_dir = "./Sensor_data"
    new_dir = f"{species}"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    parent_dir = f"./Sensor_data/{species}"
    new_dir = f"{kmer_species}"    
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
    
    output_path = f"./Sensor_data/{species}/{kmer_species}/"   
    
    Data.Generate_sequence_and_labels()
    
    # Generate data
    
    data_array,sequence = Data.Run_sensors(hws,-hws)
    
    print("Sensor data array dimensions:",data_array.shape)
    
    np.save(f"{output_path}{species}_ws{window_size}_kmers_{kmer_species}_sensor_data",data_array)
    
    np.savetxt(f"{output_path}{species}_ws{window_size}_kmers_{kmer_species}_sensor_data.csv",data_array, delimiter =",",fmt ='%s',comments='')
    
    with open(f"{output_path}{species}_ws{window_size}_kmers_{kmer_species}_sensor_sequence.txt", 'w') as f:
        f.write(sequence)
        
#%%
        
for i in range(len(kmer_species_list)):        
    
    Generate_sensor_data(species,window_size,n_sensors,kmer_species_list[i])

