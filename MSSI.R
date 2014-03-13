# libraries required for the calculation and the plotting of the MSSI
library(plyr)
library(zoo)
library(ggplot2)
library(gridExtra)
library(scales)

#read raw trajectory data, containing unique ID for each trajectory, X- and Y-coordinates and the frame
trajectories_raw <- read.csv("/Users/Frank/Documents/Postdoc/Franco_validation/Sampling_exp/JPG_videos/3 - trajectory data/ParticleLinker_Data00040.ijout.txt.txt",sep = ",")

# input the dataset
# 1) dataset containing the trajectories
# 2) specify the "unique identifier" column and the "time" column
# 3) specify the window sizes and the temporal resolution (i.e. granulosity) for which you want to calculate the SI

calculate_MSSI <- function(data,uniqueID="traj",time="frame",window_size,granulosity){

# rename columns in data according to specification
colnames(data)[colnames(data) == paste(uniqueID)] <- "uniqueID"
colnames(data)[colnames(data) == paste(time)] <- "time"

# specify resolution of the trajectory (i.e. granulosity)
temporal_simplification <- function(df,phi){
  df <- subset(df, time == min(time) | time %% phi  == 0  | time == max(time))
  return(df)
}

for (i in 1:length(granulosity)){

  # simplify trajectories according to granulosity
  trajectories <- ddply(data, .(uniqueID),function(x){temporal_simplification(x,granulosity[i])})

for (j in 1:length(window_size)){

    #make sure that all trajectories have at least as many fixes as window size 
    length <- ddply(trajectories, .(uniqueID), transform, N = length(x))
    trajectories <- trajectories[length$N > window_size[j], ]
    
    # specify rolling diff function to calculate the displacement between subsequent x or y coordinates
    roll_diff <- function(x) rollapply(x, 2, function(x) diff(x), by.column=F, fill = NA, align = "center")
    
    # run rolling diff function per trajectory
    diff_x <- ave(trajectories$x, trajectories[c("uniqueID")], FUN = roll_diff)
    diff_y <- ave(trajectories$y, trajectories[c("uniqueID")], FUN = roll_diff)
    
    # merge output with id and time
    id <- as.data.frame(cbind(trajectories$uniqueID,trajectories$time, row.names=NULL))
    names(id) <- c("uniqueID","time")
    disp <- as.data.frame(cbind(id,diff_x,diff_y))
    
    # calculate displacement based on diffs in x and y for subsequent fixes
    disp$disp <- sqrt(disp$diff_x^2+disp$diff_y^2)
    
    # use rollapply to sum displacement into gross displacement for each trajectory
    gd_extract <- function(x) rollapply(x, window_size[j], sum, fill=NA, by.column=F, align = "center")
    gd <- ave(disp$disp, disp[c("uniqueID")], FUN = gd_extract)
    #merge back with id and time
    gd <- as.data.frame(cbind(gd,id))
    
    # specify rolling diff function between first and last observation (to calculate net displacement)
    roll_diff_window <- function(x) rollapply(x, width=window_size[j]+1, function(x) diff(x,(window_size[j]))^2, fill = NA, align="center", by.column=FALSE)
    
    # apply rolling diff function for x and y-coordinates
    x_net <- ave(trajectories$x, trajectories[c("uniqueID")], FUN = roll_diff_window)
    y_net <- ave(trajectories$y, trajectories[c("uniqueID")], FUN = roll_diff_window)
    
    # calculate net displacement
    nd <- sqrt(x_net+y_net)
    #merge with id
    nd <- as.data.frame(cbind(nd,id))
    
    # shift necessary for correct association of net to gross displacement (because gd is shifted)
    if (window_size[j] %% 2 == 0) nd$time <- nd$time-granulosity[i]
    
    #merge net and gross displacement together and calculate SI
    SI <- merge(gd,nd,by=c("uniqueID","time"))
    SI <- transform(SI, SI = nd/gd)
    SI <- SI[order(SI$uniqueID,SI$time), ]
    
    # rename columns as in original dataset
    colnames(SI)[colnames(SI) == "uniqueID"] <- paste(uniqueID)
    colnames(SI)[colnames(SI) == "time"] <- paste(time)
    
    # add information on window size and granulosity to results
    SI$granulosity <- granulosity[i]
    SI$window_size <- window_size[j]

    if (i==1 & j == 1) SI_full <- SI
    if (i>=1 & j>1) SI_full <- rbind(SI_full,SI)
    }
}
return(SI_full)
}

# function call
MSSI <- calculate_MSSI(trajectories_raw,uniqueID="traj",time="frame",2:100,1:5)

# function to plot the trajectory and the corresponding MSSI
plot_MSSI <- function(raw_traj,data,uniqueID="traj",time="frame",random=T,N_traj=10,trajectory_select=select_traj){

# rename columns in data according to specification
colnames(data)[colnames(data) == paste(uniqueID)] <- "uniqueID"
colnames(data)[colnames(data) == paste(time)] <- "time"

if (random){select_traj <- sample(data$uniqueID,N_traj,replace=F)}

  for (k in 1:length(select_traj)){
  traj <-  subset(data, uniqueID == select_traj[k])

  traj_MSSI <- ggplot(traj, aes(time, window_size)) +
               geom_tile(data=traj, aes(width=as.numeric(granulosity),height=1,fill = SI)) + 
               scale_fill_gradientn(colours = c("red","cyan","black"), guide = "colourbar", limits=c(0,1))+
               theme(legend.position="bottom")+
               facet_wrap(~granulosity,ncol=4)

  traj_plot <-  subset(raw_traj, get(uniqueID) == select_traj[k])
                gg_traj <- ggplot(traj_plot, aes(x,y,label=frame))+ 
                geom_text(size=3)
  
  print(grid.arrange(gg_traj,traj_MSSI))

}
}

# call to plot function
plot_MSSI(trajectories_raw,MSSI,uniqueID="traj",time="frame",random=T,N_traj=10)

