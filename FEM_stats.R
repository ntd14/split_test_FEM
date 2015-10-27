rm(list = ls())
library('ggplot2')
df <- read.csv('sample_data_var.csv') 
library("plyr")

df$strain = df$mes*25/(1.74*200^2)*1000000
ggplot(df, aes(sample_type, strain))+ geom_boxplot()
dev.new()
ggplot(df, aes(as.factor(sample_num), strain))+ geom_boxplot()
dev.new()
hist(df$strain)
print("strain mean") #2012 =  2014 = 342.6548 2013 = 359.3972
print(mean(df$strain))
print("strain std dev") #2012 =  2014 = 332.0421  2013 = 382.004
print(sd(df$strain))
print("num of samples")
print(length(df$strain)/9)

df_real <- read.csv('../Harewood_stats/all_bosis_strain_data - Sheet1.csv')

df_real <- within(df_real, {
			treeID <- interaction(row, column, block)
			bcID <- interaction(block, coppice)
			strain <- opening*diameter/(1.74*slit^2)*1000000
		})

df_real <- df_real[!is.na(df_real$strain),] 

df_real_2012 = df_real[df_real$bcID == 'H2.0',]
df_real_2012 <- df_real_2012[order(df_real_2012$strain),]

median_strain = median(df_real_2012$strain) 
print(median_strain) # 2012 = 304.5977 
median_ind = (length(df_real_2012$strain)/2)
upper_sd_ind = round(median_ind*0.68) + median_ind

strain_sd = df_real_2012$strain[upper_sd_ind] - median_strain
print(strain_sd)	#2012 = 598.6299


df_real_2014 = df_real[df_real$bcID == 'H2.1',]
df_real_2014 <- df_real_2014[order(df_real_2014$strain),]

median_strain = median(df_real_2014$strain) 
print(median_strain) # 2014 = 942.8477
median_ind = (length(df_real_2014$strain)/2)
upper_sd_ind = round(median_ind*0.68) + median_ind

strain_sd = df_real_2014$strain[upper_sd_ind] - median_strain
print(strain_sd)	#2014 = 586.0412



