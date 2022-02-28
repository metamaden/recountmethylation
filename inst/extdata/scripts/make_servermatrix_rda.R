#!/usr/bin/env R

# Author: Sean Maden
# 
# Make the servermatrix table, to load where servermatrix() fails.

# get the table
sm <- servermatrix()
# get the ts as max ts observed
tsv <- gsub(".*_|\\..*", "", sm[,1])
tsv <- gsub(".*\\-.*", "", tsv)
tsv <- tsv[!tsv == ""]
max.ts <- max(as.numeric(tsv))
# store the new file
sm.fname <- paste0("servermatrix_",max.ts,".rda")
save(sm, file = file.path(sm.fname))