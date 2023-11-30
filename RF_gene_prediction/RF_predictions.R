### Set working directory ###

setwd("/home/.../RF_gene_prediction")

### For test species, enter the name of the target fasta file wihthout file extension

test_species <- "Drosophila_2L_minimal"

### RF parameters ###

Trained_model_name <- "10S_forest"

### Load trained RFs ###

library(randomForest)

Test_rf <- function(kmer_sp,test_sp,input_name){

load(file = paste("./Trained_models/",kmer_sp,"_",input_name,"_rf.rda", sep = "", collapse = ""))

sdata <- read.csv(paste("./Sensor_data/",test_sp,"/",kmer_sp,"/",test_sp,"_ws99_kmers_",kmer_sp,"_sensor_data.csv", sep = "", collapse = ""), header = FALSE, sep = ",")

print(model)

#use fitted model to predict class values of new data

test_data = predict(model, newdata = sdata)

write.csv(test_data, file = paste(kmer_sp,"_",test_sp,"_",input_name,"_test_results.csv", sep = "", collapse = ""))

}

kmer_species_list <- list("Drosophila_2L")

test_species_list <- list(test_species)

#test_species_list <- list("Human_18","Human_19","Human_20","Human_22")

#test_species_list <- list("Worm_II","Worm_III","Worm_IV","Worm_V")

#test_species_list <- list("Arabidopsis_2","Arabidopsis_3","Arabidopsis_4","Arabidopsis_5")

#test_species_list <- list("Saccharomyces_VII","Saccharomyces_XII","Saccharomyces_XIII","Saccharomyces_XV")

for (kmer_sp in kmer_species_list){

for (test_sp in test_species_list){

Test_rf(kmer_sp,test_sp,Trained_model_name)

}
}


