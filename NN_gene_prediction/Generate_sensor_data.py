#!/usr/bin/env python
# coding: utf-8

#%% 

import numpy as np
from collections import Counter
import random
import os

#%%

################################# PARAMETERS #################################

species_list = ["Drosophila_2L_minimal"]
chromosome_list = ["2L"]
kmer_species_list = ["Drosophila_2L_minimal"]
kmer_chromosome_list = ["2L"]
window_size = 69
n_sensors = 16
kmer_training_multiplier = 0.001

############################### END PARAMETERS ###############################

############################## Sensor functions ###############################
# All sensors are functions that take strings of letters (A,C,G,T) as input 
# and output numbers which are the features used by the model.

###################### Start codon sensor: right
# Output increases as start codon grows closer to right of target
# Maximum value attained when start codon is right next to target base

# k: distance to consider, also maximum value of sensor

def Start_codon_sensor_right(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(midpoint,len(seq)-3):
        d0 = seq[i:i+3]
        if d0.upper() == "ATG":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d

###################### Start codon sensor: left
# Output increases as start codon grows closer to left of target
# Maximum value attained when start codon is right next to target base

# k: distance to consider, also maximum value of sensor


def Start_codon_sensor_left(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(0,midpoint):
        d0 = seq[i:i+3]
        if d0.upper() == "ATG":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d


###################### Stop codon sensor: right
# Output increases as stop codon grows closer to right of target
# Maximum value attained when stop codon is directy next to target base

# k: distance to consider, also maximum value of sensor

def Stop_codon_sensor_right(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(midpoint,len(seq)-3):
        d0 = seq[i:i+3]
        if d0.upper() == "TAA" or d0.upper() == "TAG" or d0.upper() == "TGA":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d

##################### Stop codon sensor: left
# Output increases as stop codon grows closer to left of target
# Maximum value attained when stop codon is directy next to target base

# k: distance to consider, also maximum value of sensor

def Stop_codon_sensor_left(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(0,midpoint):
        d0 = seq[i:i+3]
        if d0.upper() == "TAA" or d0.upper() == "TAG" or d0.upper() == "TGA":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d

##################3# Donor splice site sensor: right
# Output increases as donor site grows closer to right of target
# Maximum value attained when splice site is directy next to target base

# k: distance to consider, also maximum value of sensor

def Donor_splice_site_sensor_right(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(midpoint,len(seq)-3):
        d0 = seq[i:i+3]
        if d0.upper() == "GT":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d

####################### Donor splice site sensor: left
# Output increases as donor site grows closer to left of target
# Maximum value attained when splice site is directy next to target base

# k: distance to consider, also maximum value of sensor

def Donor_splice_site_sensor_left(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(0,midpoint):
        d0 = seq[i:i+3]
        if d0.upper() == "GT":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d

################ Acceptor splice site sensor: right
# Output increases as accepor site grows closer to right of target
# Maximum value attained when splice site is directy next to target base

# k: distance to consider, also maximum value of sensor

def Acceptor_splice_site_sensor_right(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(midpoint,len(seq)-3):
        d0 = seq[i:i+3]
        if d0.upper() == "AG":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d

################### Acceptor splice site sensor: left
# Output increases as accepor site grows closer to left of target
# Maximum value attained when splice site is directy next to target base

# k: distance to consider, also maximum value of sensor

def Acceptor_splice_site_sensor_left(seq):
    midpoint = int((len(seq) - 1)/2) # center of window
    d = 0
    for i in range(0,midpoint):
        d0 = seq[i:i+3]
        if d0.upper() == "AG":
            if (midpoint-abs(midpoint-i)) > d: 
                d = midpoint-abs(midpoint-i)
    return d

############### Voss encoder and DFT period 3 sensor
# outputs period 3 signal from input sequence
# w0: input sequence

def Voss_encoder(seq):
    uA = []
    uT = []
    uC = []
    uG = []
    for i in seq:
        if i.upper() == 'A':
            uA.append(1)
            uT.append(0)
            uC.append(0)
            uG.append(0)
        if i.upper() == 'T':
            uA.append(0)
            uT.append(1)
            uC.append(0)
            uG.append(0)
        if i.upper() == 'C':
            uA.append(0)
            uT.append(0)
            uC.append(1)
            uG.append(0)
        if i.upper() == 'G':
            uA.append(0)
            uT.append(0)
            uC.append(0)
            uG.append(1)
        else:
            uA.append(0)
            uT.append(0)
            uC.append(0)
            uG.append(0)
    return uA,uT,uC,uG

def DFT_L3_sensor(seq):
    wA,wT,wC,wG = Voss_encoder(seq)
    i = np.sqrt(-1+0j)
    Pi = 3.14159
    N = len(wA)  # N = length of window
    k = N/3
    sumA = 0
    sumT = 0
    sumC = 0
    sumG = 0
    for n in range(0,N): # n = window position index
        exp0 = np.exp(-i*2*Pi*(k)*(n)/N)
        sumA = sumA + wA[n]*exp0
        sumT = sumT + wT[n]*exp0
        sumC = sumC + wC[n]*exp0
        sumG = sumG + wG[n]*exp0
        
    sum0 = ((sumA*np.conjugate(sumA)) + (sumT*np.conjugate(sumT)) + 
            (sumC*np.conjugate(sumC)) + (sumG*np.conjugate(sumG)))
    
    return sum0.real

################# GC content sensor
# seq: input sequence in string format
# outputs GC content of input sequence

def GC_content_sensor(seq):
    gc_count = 0
    for i in seq:
        if i.upper() == 'C':
            gc_count +=1
        if i.upper() == 'G':
            gc_count +=1
    return gc_count/len(seq)


################# Kmer sensor
    
class Kmer_sensor():
    
    def __init__(self, coding_freqs=None, noncoding_freqs=None,k=4):
        
        self.coding_freqs = Counter()
        self.noncoding_freqs = Counter()
        self.tuple_size = k

    def Train_kmer_sensor(self,seq,dnalabels,ws0,kmer_training_multiplier):

        hws = int((ws0-1) / 2) # half window size. Set to zero to ignore cds/noncds boundaries.
        
#        dnaseq = seq

        feature_list = []

        non_feature_list = []

        for i in range(hws,len(dnalabels)-hws):
            if dnalabels[i] == 1:
                if 0 not in dnalabels[i-hws:i+hws]:
                    feature_list.append(i)
            if dnalabels[i] == 0:
                if 1 not in dnalabels[i-hws:i+hws]:
                    non_feature_list.append(i)
        
        training_number = int( kmer_training_multiplier * (4**self.tuple_size) * 1000 * self.tuple_size * 2 / ws0 )
        
        x_train0 = []
        y_train0 = []

        for i in range(0,training_number,2):
            target = random.choice(feature_list)
            x_train0.append(seq[(target-hws):(target+hws+1)].upper())
            y_train0.append(1)

            target = random.choice(non_feature_list)
            x_train0.append(seq[(target-hws):(target+hws+1)].upper())
            y_train0.append(0)

        # For each row in data
        for row, label in zip(x_train0, y_train0):

            # Stringify
            sub_seq = "" 
            sub_seq = sub_seq.join(row)

            if label == 0:
                self.noncoding_freqs += self._count_kmers(sub_seq, self.tuple_size)
            elif label == 1:
                self.coding_freqs += self._count_kmers(sub_seq, self.tuple_size)
            else:
                raise ValueError('Labels must be 0 or 1')        
                
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
    
###################### Shannon entropy sensor

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

###################### Sensor sweep code
# Outputs a vector of sensor values
# cdi = sweep dna start coordinate
# cdf = sweep dna stop coordinate
# seq: input dna sequence from fasta
# ws: window size

def Sensor_sweep(seq, target_nt,ws,k6_sensor,k5_sensor,k4_sensor,k3_sensor,k2_sensor,nsensors):
    output_array = np.zeros(nsensors)
    r = int((ws-1)/2) # half window size
    row = 0
    i = target_nt
    target_sequence = seq[i-r:i+r+1]
    output_array[row] = DFT_L3_sensor(target_sequence); row = row + 1;
    output_array[row] = GC_content_sensor(target_sequence); row = row + 1;
    output_array[row] = k6_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k5_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k4_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k3_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = k2_sensor.encode(target_sequence); row = row + 1;
    output_array[row] = Donor_splice_site_sensor_right(target_sequence); row = row + 1;
    output_array[row] = Acceptor_splice_site_sensor_right(target_sequence); row = row + 1;
    output_array[row] = Start_codon_sensor_right(target_sequence); row = row + 1;
    output_array[row] = Stop_codon_sensor_right(target_sequence); row = row + 1;
    output_array[row] = Donor_splice_site_sensor_left(target_sequence); row = row + 1;
    output_array[row] = Acceptor_splice_site_sensor_left(target_sequence); row = row + 1;
    output_array[row] = Start_codon_sensor_left(target_sequence); row = row + 1;
    output_array[row] = Stop_codon_sensor_left(target_sequence); row = row + 1;
    output_array[row] = Entropy_sensor(target_sequence); row = row + 1;
    return output_array

#%%
    
############################## Data generator ##############################

# fasta file and gff file must be in the same directory as this code
# Default values of class variables are set to use a small test dataset
    
class Data_generator():
    
    def __init__(self,ws,ns,ktm):
        
        self.fasta_input = "Drosophila_2L_minimal.fna" # my filename for the input fasta file
        self.gff_input = "Drosophila_2L_minimal.gff3" # my filename for the input gff file
        self.kmer_fasta_input = "Drosophila_2L_minimal.fna" # my filename for the input fasta file
        self.kmer_gff_input = "Drosophila_2L_minimal.gff3" # my filename for the input gff file
        
        self.feature = "CDS" # sets type of feature to consider for positive samples
        self.nsensors = ns
        self.chromosome = "2L"
        self.kmer_chromosome = "2L"
        self.window_size = ws
        self.kmer_training_multiplier = ktm
        
# Reads a gff file and stores the start and stop coordinates for each feature of the given type
# Reads a fasta file to generate a DNA string object for use by other functions
# Generates a label string the same length as the DNA string object
# Label string has '1' where there is a feature of the given type (eg. CDS)
# Label string has '0' where there is no feature of the given type (eg. CDS)
# Locations of features (eg. CDS) are extracted from inout gff file
        
    def Generate_sequence_and_labels(self):

        self.slice0 = []
        
        with open(self.gff_input) as gff_file:
            gff = gff_file.readlines()
            for line in gff:
                if line.split()[0] == self.chromosome:
                    if line.split()[2] == self.feature:
                        if line.split()[6] == "+":
                            self.slice0.append((line.split()[3],line.split()[4]))
                                
                            genome = open(self.fasta_input)
                    
                            lines = genome.readlines()
                            DNASeq0 = []
                            for i in range (0, len(lines)):
                                if lines[i][0:1] != ">":
                                    DNASeq0.append(lines[i].strip("\n"))
                            genome.close()
                            
                            self.DNASeq = ''.join(DNASeq0)
                            
                            self.DNALabels = np.zeros((len(self.DNASeq)))
                    
                            for j in range(0,len(self.slice0)):
                                start = int(self.slice0[j][0])
                                stop = int(self.slice0[j][1])
                                if stop < len(self.DNASeq):
                                    for k in range(start,stop):
                                        self.DNALabels[k] = 1

#%%
                                 
    def Generate_kmer_sequence_and_labels(self):

        self.kmer_slice0 = []
        
        with open(self.kmer_gff_input) as gff_file:
            gff = gff_file.readlines()
            for line in gff:
                if line.split()[0] == self.kmer_chromosome:
                    if line.split()[2] == self.feature:
                        if line.split()[6] == "+":
                            self.kmer_slice0.append((line.split()[3],line.split()[4]))
                                
                            genome = open(self.kmer_fasta_input)
                    
                            lines = genome.readlines()
                            DNASeq0 = []
                            for i in range (0, len(lines)):
                                if lines[i][0:1] != ">":
                                    DNASeq0.append(lines[i].strip("\n"))
                            genome.close()
                            
                            self.kmer_DNASeq = ''.join(DNASeq0)
                            
                            self.kmer_DNALabels = np.zeros((len(self.kmer_DNASeq)))
                    
                            for j in range(0,len(self.kmer_slice0)):
                                start = int(self.kmer_slice0[j][0])
                                stop = int(self.kmer_slice0[j][1])
                                if stop < len(self.kmer_DNASeq):
                                    for k in range(start,stop):
                                        self.kmer_DNALabels[k] = 1
                                           
                    
####### Function to train kmer sensors

    def Train_kmer_sensors(self):

        self.kmer_6_sensor = Kmer_sensor(k=6)
        self.kmer_6_sensor.Train_kmer_sensor(self.kmer_DNASeq,self.kmer_DNALabels,self.window_size,
                                             self.kmer_training_multiplier)
        
        self.kmer_5_sensor = Kmer_sensor(k=5)
        self.kmer_5_sensor.Train_kmer_sensor(self.kmer_DNASeq,self.kmer_DNALabels,self.window_size,
                                             self.kmer_training_multiplier)
        
        self.kmer_4_sensor = Kmer_sensor(k=4)
        self.kmer_4_sensor.Train_kmer_sensor(self.kmer_DNASeq,self.kmer_DNALabels,self.window_size,
                                             self.kmer_training_multiplier)
        
        self.kmer_3_sensor = Kmer_sensor(k=3)
        self.kmer_3_sensor.Train_kmer_sensor(self.kmer_DNASeq,self.kmer_DNALabels,self.window_size,
                                             self.kmer_training_multiplier)
        
        self.kmer_2_sensor = Kmer_sensor(k=2)
        self.kmer_2_sensor.Train_kmer_sensor(self.kmer_DNASeq,self.kmer_DNALabels,self.window_size,
                                             self.kmer_training_multiplier)

######### Run_sensors() generates sensor data with label column and nt sequence
    
    def Run_sensors(self,startnt,stopnt):
        
        labels0 = self.DNALabels[startnt:stopnt]
        
        test_seq = self.DNASeq[startnt:stopnt]
        
        test_number = len(labels0)
                
        x_test0 = np.zeros((test_number,self.nsensors))

        for i in range(0,test_number):
            target = startnt + i
            x_test0[i] = Sensor_sweep(self.DNASeq, target,self.window_size,
                                      self.kmer_6_sensor,self.kmer_5_sensor,self.kmer_4_sensor,
                                      self.kmer_3_sensor,self.kmer_2_sensor,self.nsensors)
            
        label_column = labels0.reshape(-1,1)
            
        sensor_data = np.hstack((x_test0, label_column))
                
        return sensor_data,test_seq      
    
#%%

def Generate_sensor_data(species,window_size,
                             nsensors,kmer_species,chromosome,kmer_chromosome,
                             kmer_training_multiplier):
    
    hws = int( (window_size - 1 ) / 2 )

    Data = Data_generator(window_size,nsensors,kmer_training_multiplier)
    
    fasta_path = "./fasta_files/"
    gff_path = "./gff_files/"
    
    parent_dir = "./"
    new_dir = "Sensor_data"
    path = os.path.join(parent_dir, new_dir)
    if not os.path.isdir(path):
        os.mkdir(path)
        
    Data.fasta_input = fasta_path + species + ".fna"
    Data.gff_input = gff_path + species + ".gff3"
    Data.chromosome = chromosome
    
    Data.kmer_fasta_input = fasta_path + kmer_species + ".fna"
    Data.kmer_gff_input = gff_path + kmer_species + ".gff3"
    Data.kmer_chromosome = kmer_chromosome
    
    Data.Generate_kmer_sequence_and_labels()
    Data.Train_kmer_sensors()
    print("Kmer sensor training secies:",kmer_species)
  
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
    
    data_array,sequence = Data.Run_sensors(hws,-hws)
    
    print("Sensor data array dimensions:",data_array.shape)
    
    np.save(f"{output_path}{species}_ws{window_size}_kmers_{kmer_species}_sensor_data",data_array)
    
    np.savetxt(f"{output_path}{species}_ws{window_size}_kmers_{kmer_species}_sensor_data.csv",data_array, delimiter =",",fmt ='%s',comments='')
    
    with open(f"{output_path}{species}_ws{window_size}_kmers_{kmer_species}_sensor_sequence.txt", 'w') as f:
        f.write(sequence)
        
#%%
        
for i in range(len(species_list)):
    for j in range(len(kmer_species_list)):        
        
        Generate_sensor_data(species_list[i],window_size,n_sensors,
                             kmer_species_list[j],chromosome_list[i],
                             kmer_chromosome_list[j],kmer_training_multiplier)

