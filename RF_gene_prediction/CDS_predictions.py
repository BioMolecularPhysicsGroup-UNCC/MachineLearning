#!/usr/bin/env python
# coding: utf-8

#%% 

import numpy as np
import csv

#%%

################################# USER PARAMETERS #################################

### For test species, enter the name of the target fasta file wihthout file extension

test_species = "Drosophila_2L_minimal"

### For path, enter working directory

path = "/home/.../RF_gene_prediction/"

################################# OTHER PARAMETERS #################################

max_gt = 10
max_stop = 2
min_cds_length = 40
initial_cds_threshold = 0.9
max_cds_threshold = 1.0
threshold_step = 0.1

kmer_species_list = ["Drosophila_2L"]

model_name = "10S_forest"

#%%

def Generate_prediction_array(kmer_species_list,test_species):
    
    prediction_vectors = []
    
    for train_sp in kmer_species_list:        
    
        pred_path = train_sp + "_" + test_species + "_" + model_name + "_test_results.csv"
        
        with open(pred_path,'r') as p:
            pred_list = [row[1] for row in csv.reader(p, delimiter=',')]
            pred = np.zeros(len(pred_list)-1)
            for i in range(0,len(pred_list)-1):
                pred[i] = int(pred_list[i+1])
                
            prediction_vectors.append(pred)
            
    n_predictions = prediction_vectors[0].shape[0]
            
    combined_predictions = np.zeros(n_predictions)
    
    num_votes = len(prediction_vectors) # number of votes per sample
    
    for i in range(0,n_predictions):
        
        sum_vote = 0
        
        for j in range(0,num_votes):            
            sum_vote += prediction_vectors[j][i]
       
        if sum_vote/num_votes >= 0.5:
            combined_predictions[i] = 1
            
    return combined_predictions
            

#pred_array = Generate_prediction_array(training_species_list,"Drosophila_2L_minimal")

#%%


def Adjust_predictions(kmer_species_list,test_gene,seq):

    with open(seq,'r') as f:
        sequence = f.read()

    seq_len = len(sequence)

    pred = Generate_prediction_array(kmer_species_list,test_gene)

    ### nucleeotide level results

    total_cds = 0

    for i in range(0,seq_len-1):
        if pred[i] == 1:
            total_cds += 1

    base_cds_fraction = total_cds/seq_len

    ### scan exon sites

    ### maximum possible number of stop codons in single exon which are not in_frame = 3

    potential_exon_sites = []

    for i in range(0,seq_len-3):
        if sequence[i:i+3].upper() == "ATG":
            stop_count = 0
            for j in range(i+min_cds_length,seq_len-3):
                if sequence[j:j+3].upper() in ["TAG","TAA","TGA"]:
                    if stop_count > max_stop:
                        break
                    stop_count += 1
                    if len(sequence[i:j+3]) % 3 == 0:
                        cds = sum(pred[i:j])
                        if cds/(j-i) > initial_cds_threshold:
                            potential_exon_sites.append((i,j,cds/(j-i),"Single"))
                        break ### break after finding a single in_frame stop codon

            gt_count = 0
            for j in range(i+min_cds_length,seq_len-3):
                if sequence[j:j+3].upper() in ["TAG","TAA","TGA"]:
                    if len(sequence[i:j+3]) % 3 == 0: ### check for one in_frame stop codon
                        break                
                if sequence[j:j+2].upper() == "GT":
                    if gt_count > max_gt:
                        break
                    cds = sum(pred[i:j])
                    gt_count += 1
                    if cds/(j-i) > initial_cds_threshold:
                        potential_exon_sites.append((i,j,cds/(j-i),"Start"))

        if sequence[i:i+3].upper() in ["AAG","CAG","TAG"]:
            stop_count = 0
            stop_codon_sites = []
            for j in range(i+min_cds_length,seq_len-3):
                in_frame = 0
                if sequence[j:j+3].upper() in ["TAG","TAA","TGA"]:
                    if stop_count > max_stop:
                        break
                    stop_count += 1
                    if len(stop_codon_sites) > 0:
                        for k in stop_codon_sites:
                            if (j - k) % 3 == 0: ### check for in frame stop codons
                                in_frame = 1
                    stop_codon_sites.append(j)
                    if in_frame == 0:
                        cds = sum(pred[i:j])
                        if cds/(j-i) > initial_cds_threshold:
                            potential_exon_sites.append((i,j,cds/(j-i),"Stop"))

            gt_count = 0
            stop_codon_sites = []
            stop_frames = 0
            for j in range(i+min_cds_length,seq_len-3):
                in_frame = 0
                if sequence[j:j+3].upper() in ["TAG","TAA","TGA"]:
                    if len(stop_codon_sites) > 0:
                        for k in stop_codon_sites:
                            if (j - k) % 3 == 0: ### check for in frame stop codons
                                in_frame = 1
                    if in_frame == 0:
                        stop_frames += 1
                    if stop_frames > 2: ### maximum 3 out of frame stop codons possible
                        break
                    stop_codon_sites.append(j)
                if sequence[j:j+2].upper() == "GT":
                    if gt_count > max_gt:
                        break
                    cds = sum(pred[i:j])
                    gt_count += 1
                    if cds/(j-i) > initial_cds_threshold:
                        potential_exon_sites.append((i,j,cds/(j-i),"Internal"))


    ### count potential exon sites backwards from largest to smallest
    ### Skip over small potential exons if larger predicted exon encompases it

    for threshold in np.arange(initial_cds_threshold,max_cds_threshold,threshold_step):

        def filter_sites(potential_exon_site):
            internal_cds_fraction = potential_exon_site[2]
            cds_ratio = internal_cds_fraction / base_cds_fraction
            if cds_ratio > threshold:
                return potential_exon_site

        f_potential_exon_sites = list(filter(filter_sites,potential_exon_sites))

        predictions = np.zeros((seq_len))

        exon_lengths = []
        intron_lengths = []

        single_exons = 0
        start_exons = 0
        stop_exons = 0
        internal_exons = 0

        total_cds = 0

        print("--------------------------------------------------")
        print("Max CDS length:",max_gt,"GT sites",max_stop,"stop sites")
        print("CDS fraction threshold:",np.round(threshold,1))

        predicted_exon_sites = []

        exon_start = 0
        best_exon = (0,0,threshold,"null")
        best_exon_length = 0

        for k in range(0,len(f_potential_exon_sites)):

            exon_start = f_potential_exon_sites[k][0]
            exon_stop = f_potential_exon_sites[k][1]
            exon_length = exon_stop - exon_start

            if exon_start > int(best_exon[0]):
                if best_exon[0] != 0:
                    predicted_exon_sites.append(best_exon)
#                    print(best_exon)
                    for m in range(int(best_exon[0]),int(best_exon[1])):
                        predictions[m] = 1
                    exon_lengths.append(best_exon_length)
                    if len(predicted_exon_sites) > 1:
                        intron_lengths.append(int(best_exon[0]) - predicted_exon_sites[-2][1])
                    best_exon = (0,0,threshold,"null") # reset for new CDS group
                    best_exon_length = 0


            if exon_length > best_exon_length:
                best_exon = f_potential_exon_sites[k]
                best_exon_length = exon_length

### Add last predicted exon

        if best_exon[0] != 0:
            predicted_exon_sites.append(best_exon)
#            print(best_exon)
            for m in range(int(best_exon[0]),int(best_exon[1])):
                predictions[m] = 1
            exon_lengths.append(best_exon_length)
            if len(predicted_exon_sites) > 1:
                intron_lengths.append(int(best_exon[0]) - predicted_exon_sites[-2][1])

        ### count predicted exon types

        for p in range(0,len(predicted_exon_sites)):
            if predicted_exon_sites[p][3] == "Single":
                single_exons += 1
            if predicted_exon_sites[p][3] == "Start":
                start_exons += 1
            if predicted_exon_sites[p][3] == "Stop":
                stop_exons += 1
            if predicted_exon_sites[p][3] == "Internal":
                internal_exons += 1
                
        for i in range(len(predicted_exon_sites)):
            
            print(test_species,"CDS",predicted_exon_sites[i][0],predicted_exon_sites[i][1],"+")
                
        
        mean_exon_length = 0

        if len(exon_lengths) > 0:
            mean_exon_length = np.mean(exon_lengths[:])
            
        print("Mean CDS length:", mean_exon_length)

        mean_intron_length = 0

        if len(intron_lengths) > 0:
            mean_intron_length = np.mean(intron_lengths[:])
            
        print("Mean intron length:",mean_intron_length)
        
        print(len(predicted_exon_sites),"CDS regions predicted")

#%%

def GPU():
        
    print("TEST SPECIES RESULTS:",test_species)
    
    print("BEGIN TRAINING SPECIES LIST")
    
    for species in kmer_species_list:
        print(species)
            
    print("END TRAINING SPECIES LIST")

    seq_path = path + "Sensor_data/" + test_species + "/" + kmer_species_list[0] + "/" + test_species + "_ws99_kmers_" + kmer_species_list[0] + "_sensor_sequence.txt"    
    Adjust_predictions(kmer_species_list,test_species,seq_path)
    
    print("END TEST SPECIES RESULTS:",test_species)
    
    
#%%
    
GPU()
