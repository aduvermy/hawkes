rm(list=ls())

library(fs)
library(stringr)
library(tidyverse)

##HAWKES PARAMS
exec="~/hawkes/hawkes"   #path to hawkes executable
delta  = 1000 ## taille des fenêtres sur lesquelles on regarde les interactions
K      = 20  ## nombre de fenêtres
kernel = "heterogeneous_interval" ## "none" = cas ponctuel, "homogeneous_interval" intervalles de longueur fixée à la longueur moyenne, "heterogeneous_interval" intervalles de longueurs variables (+ long)
lambda = 1 ## paramètre de pénalisation, plus il est petit moins on met de constantes à 0


##BED INPUT
path2BED<-"BED_directory"


my_bed_files <- dir_ls(path = path2BED ,glob = "*.bed")
my_bed_files

##OUTPUT DIRECTORY (without last /)

outdir="output directory"

dir_create(outdir)

#BUILD OUTPUT SIMPLE NAME
pattern <- ".+/(.+)\\.\\w+$"
simpleName_output<- sub(pattern, replacement = "\\1", my_bed_files)%>%
                        str_c(collapse = "_")%>%
                        str_c("K",K,"delta",delta,"kernel",kernel,'lambda',lambda,sep = "_")

source("hawkes_function.R") ##path to hawkes functions
simpleName_output

##RUN HAWKES
forward=hawkes( "f" , my_bed_files, simpleName_output)
backward=hawkes( "b" , my_bed_files, simpleName_output)
f=getparam(forward)
b=getparam(backward)
dta=merge_fb(f,b)


##RESULTS
file_out<-  str_c(outdir,simpleName_output,sep ="/")%>%
                str_c('svg',sep=".")
file_out
svg(filename=file_out, width=30, height=30)
    plot_res(my_bed_files,dta)
dev.off()
