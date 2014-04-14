# libraries required for the calculation and the plotting of the MSSI
library(dplyr)
library(zoo)
library(ggplot2)
library(gridExtra)
#library(scales)
library(tcltk)

#read raw trajectory data, containing unique ID for each trajectory, X- and Y-coordinates and the frame
trajectory.data.full <- read.table("/Users/Frank/Documents/Postdoc/Franco_validation/Sampling_exp/PNG_videos/5 - merged data/MasterData.csv",header=T,row.names=NULL,sep=",",stringsAsFactor=F)
trajectory.data <- trajectory.data.full[trajectory.data.full$file == "data00001", ]
trajectory.data <- trajectory.data[order(trajectory.data$file,trajectory.data$trajectory,trajectory.data$frame), ]
#create unique ID consisting of trajectory ID and file
traj <- paste(trajectory.data$file,trajectory.data$trajectory,sep="-")
trajectory.data <- cbind(trajectory.data,traj)

#' Function to calculate the Multiscale Straightness Index
#' 
#' Takes a dataframe with the X and Y coordinates of a trajectory and calculates the MSSI for different window sizes and 
#' granulosities
#' @param data A dataframe containing the X and Y coordinates
#' @param uniqueID The unique identifier for each trajectory in the dataframe
#' @param time The column containing the information on the time for each fix
#' @param window_size A numeric value or vector specifying the window sizes over which the MSSI should be calculated
#' @param granulosity A numeric value or vector specifying the temporal resolution of the data before the MSSI is calculated
#' @return A dataframe containing the MSSI for each fix in addition to the input data
#' @import dplyr zoo tcltk
#' @export

calculate_MSSI <- function(data,uniqueID="traj",time="frame",window_size,granulosity){

original_id <- eval(get(uniqueID))
  
# rename columns in data according to specification
colnames(data)[colnames(data) == paste(uniqueID)] <- "uniqueID"
colnames(data)[colnames(data) == paste(time)] <- "time"
colnames(data)[colnames(data) == "X"] <- "x"
colnames(data)[colnames(data) == "Y"] <- "y"

data <- data[,c("uniqueID","time","x","y")]

# part to account for character values as identifiers
data$uniqueID <- as.numeric(as.factor(data$uniqueID))
original_id <- as.data.frame(cbind(original_id,data$uniqueID))
original_id <- unique(original_id)
rownames(original_id) <- NULL
colnames(original_id) <- c("original_id","uniqueID")
rownames(data) <- NULL

# create progress bar
pb <- tkProgressBar(title = "progress bar", min = 0, max = length(window_size), width = 300)

for (i in 1:length(granulosity)){

# simplify trajectories according to granulosity
trajectories <- as.data.frame(data %.%
                group_by(uniqueID) %.%
                filter(time == min(time) | time %% granulosity[i]  == 0  | time == max(time))%.%
                arrange(uniqueID,time))

for (j in 1:length(window_size)){
  
   Sys.sleep(0.1)
   setTkProgressBar(pb, j, label=paste(round(j/length(window_size)*100, 0), "% done of run ",i))
  
    #make sure that all trajectories have at least as many fixes as window size 
    length <- as.data.frame(trajectories %.%
              group_by(uniqueID) %.%
              mutate(N = length(x)))
   
    trajectories <- trajectories[length$N > window_size[j], ]

    if(nrow(trajectories)==0){ next
    
    } else {
    
    # specify rolling diff function to calculate the displacement between subsequent x or y coordinates
    roll_diff <- function(x) rollapply(x, 2, function(x) diff(x), by.column=F, fill = NA, align = "center")
 
    # run rolling diff function per trajectory
    disp <- as.data.frame(trajectories %.%
                         group_by("uniqueID") %.%
                         transform(diff_x = roll_diff(x), diff_y=roll_diff(y)))
  
    # calculate displacement based on diffs in x and y for subsequent fixes
    disp$disp <- sqrt(disp$diff_x^2+disp$diff_y^2)
    disp <- disp[, c(1,2,7)] 
    
    # use rollapply to sum displacement into gross displacement for each trajectory
    gd_extract <- function(x) rollapply(x, window_size[j], sum, fill=NA, by.column=F, align = "center")

    gd <- as.data.frame(disp %.%
          group_by(uniqueID) %.%
          mutate(gd = gd_extract(disp)))

    gd$disp <- NULL

    # specify rolling diff function between first and last observation (to calculate net displacement)
    roll_diff_window <- function(x) rollapply(x, width=window_size[j]+1, function(x) diff(x,(window_size[j]))^2, fill = NA, align="center", by.column=FALSE)
    
    nd <- as.data.frame(trajectories %.%
          group_by("uniqueID") %.%
          transform(x_net = roll_diff_window(x),
                    y_net = roll_diff_window(y)))

    # calculate net displacement
    nd$nd <- sqrt(nd$x_net+nd$y_net)
    nd <- nd[, c(1,2,7)] 

    # shift necessary for correct association of net to gross displacement (because gd is shifted)
    if (window_size[j] %% 2 == 0) nd$time <- nd$time-granulosity[i]
    
    #merge net and gross displacement together and calculate SI
    SI <- merge(gd,nd,by=c("uniqueID","time"))
    SI <- transform(SI, SI = nd/gd)
    SI <- SI[order(SI$uniqueID,SI$time), ]
    SI$nd <- SI$gd <- NULL
    
    SI <- merge(SI,original_id,c("uniqueID"))
    SI$uniqueID <- NULL
    
    # rename columns as in original dataset
    colnames(SI)[colnames(SI) == "original_id"] <- paste(uniqueID)
    colnames(SI)[colnames(SI) == "time"] <- paste(time)
    
    # add information on window size and granulosity to results
    SI$granulosity <- granulosity[i]
    SI$window_size <- window_size[j]

    if (i==1 & j == 1) SI_full <- SI else SI_full <- rbind(SI_full,SI)
        }
    }
}
close(pb)
rownames(SI_full) <- NULL
return(SI_full)
}

# calculate_MSSI function call
MSSI <- calculate_MSSI(trajectory.data,uniqueID="traj",time="frame",2:3,1)

#' Function to visualize the MSSI along a trajectory 
#' 
#' Takes the dataframe returned by the calculate_MSSI function and the raw trajectory data to construct the plot
#' @param raw_traj A dataframe containing the X and Y coordinates of the trajectories on which the MSSI was calculated
#' @param data The dataframe returned by the calculate MSSI_function containing the MSSI index for each fix and different window sizes 
#' and granulosities
#' @param uniqueID Column name of the unique identifier for each trajectory in the dataframe
#' @param time Column name containing the information on the time for each fix
#' @param random Logical; whether the trajectories to be plotted are randomly sampled from the data 
#' @param N_traj A numeric value specifying the number of trajectories to be plotted
#' @return A plot showing the MSSI for each fix for different window sizes
#' @param trajectory_select Vector containing the unique IDs of the trajectories to be plotted'
#' @import ggplot2 gridExtra 
#' @export

plot_MSSI <- function(raw_traj,data,uniqueID="traj",time="frame",random=T,N_traj=10,trajectory_select=select_traj){

#rename x and y variable to make function cap-insensitive
colnames(raw_traj)[tolower(colnames(raw_traj)) == "x"] <- "x"
colnames(raw_traj)[tolower(colnames(raw_traj)) == "y"] <- "y"

# rename columns in data according to specification
colnames(data)[colnames(data) == paste(uniqueID)] <- "uniqueID"
colnames(data)[colnames(data) == paste(time)] <- "time"

colnames(raw_traj)[colnames(raw_traj) == paste(uniqueID)] <- "uniqueID"
colnames(raw_traj)[colnames(raw_traj) == paste(time)] <- "time"

if (random){select_traj <- sample(data$uniqueID,N_traj,replace=F)}

  for (k in 1:length(select_traj)){
  traj <-  subset(data, uniqueID == select_traj[k])

  traj_MSSI <- ggplot(traj, aes(time, window_size)) +
               geom_tile(data=traj, aes(width=as.numeric(granulosity),height=as.numeric(2),fill = SI)) + 
               scale_fill_gradientn(colours = c("red","cyan","black"), guide = "colourbar", limits=c(0,1))+
               theme(legend.position="bottom")+
               facet_wrap(~granulosity,ncol=4)

  traj_plot <-  subset(raw_traj, uniqueID == select_traj[k])
                gg_traj <- ggplot(traj_plot, aes(x,y,label=time))+ 
                geom_text(size=3)
  
  print(grid.arrange(gg_traj,traj_MSSI))

}
}

# call to plot function
plot_MSSI(trajectory.data,MSSI,uniqueID="traj",time="frame",random=T,N_traj=10)






