#library(devtools)
#install_github("MSSI" ,"pennekampster")
#library(MSSI)
# getwd()
#read raw trajectory data, containing unique ID for each trajectory, X- and Y-coordinates and the frame
# trajectory.data.full <- read.table("/Users/Frank/Documents/Postdoc/Pairwise_species_compare/5 - merged data/MasterData.csv",header=T,row.names=NULL,sep=",",stringsAsFactor=F)
# trajectory.data <- trajectory.data.full[trajectory.data.full$file == "data00026", ]
# trajectory.data <- trajectory.data[order(trajectory.data$file,trajectory.data$trajectory,trajectory.data$frame), ]
# #create unique ID consisting of trajectory ID and file
# traj <- paste(trajectory.data$file,trajectory.data$trajectory,sep="-")
# trajectory.data <- cbind(trajectory.data,traj)
# trajectory.data <- subset(trajectory.data, traj != "Data00026-NA")
# example_data <- subset(trajectory.data, select=c(traj,frame,X,Y))
# save(example_data,file="example_data.rda")

data(example_data)


# calculate_MSSI function call
MSSI <- calculate_MSSI(example_data,uniqueID="traj",time="frame",20,20)

# call to plot function
plot_MSSI(example_data,MSSI,uniqueID="traj",time="frame",granulosity_choosen=10,random=T,N_traj=1)



