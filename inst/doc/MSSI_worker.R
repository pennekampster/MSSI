#library(devtools)
#install_github("MSSI" ,"pennekampster")
#library(MSSI)
# getwd()
#read raw trajectory data, containing unique ID for each trajectory, X- and Y-coordinates and the frame
# library(data.table)
# trajectory.data.full <- fread("/Users/Frank/Documents/Postdoc/Pairwise_species_compare/5 - merged data/MasterData.csv")
#  trajectory.data <- trajectory.data.full[trajectory.data.full$file == "data00008", ]
#  trajectory.data <- trajectory.data[order(trajectory.data$file,trajectory.data$trajectory,trajectory.data$frame), ]
# # #create unique ID consisting of trajectory ID and file
#  traj <- paste(trajectory.data$file,trajectory.data$trajectory,sep="-")
#  trajectory.data <- cbind(trajectory.data,traj)
#  trajectory.data <- subset(trajectory.data, traj != "data00008-NA")
#  


example_data <- as.data.frame(subset(trajectory.data2, select=c(traj,frame,X,Y)))
example_data <- example_data[,c("traj","frame","X","Y")]
save(example_data,file="example_data.rda")

data(example_data)

str(example_data)
rownames(example_data) <- NULL

uniqueID

trajectory.data2 <- as.data.frame(trajectory.data)

# calculate_MSSI function call
MSSI <- calculate_MSSI(example_data,uniqueID="traj",time="frame",2:5,1)

plot(trajectory.data2$X,trajectory.data2$Y)

# call to plot function
plot_MSSI(example_data,MSSI,uniqueID="traj",time="frame",granulosity_choosen=1,random=T,N_traj=2)

example_data$traj <- as.character(example_data$traj)


ex <- trajectory.data %.%
      group_by(traj) %.%
      summarise(length = length(frame))%.%
      filter(length>400)

trajectory.data2 <- trajectory.data[trajectory.data$traj %in% ex$traj, ]


