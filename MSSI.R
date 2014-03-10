# libraries required for the calculation and the plotting of the MSSI
library(plyr)
library(zoo)
library(ggplot2)
library(gridExtra)
library(scales)

#read raw trajectory data, containing unique ID for each trajectory, X- and Y-coordinates and the frame
trajectories_raw <- read.csv("/Users/Frank/Desktop/video_compress/data/test_videos/3 - trajectory data/ParticleLinker_Data00072.ijout.txt.txt",sep = ",")
#trajectories <- read.csv("/Users/Frank/Dropbox/test_traj.csv",sep = ";")

# specify window size (NB they must be > 3) and granulosity 
granulosity <- seq(1,20,2)
window_size <- 2:20

for (i in 1:length(granulosity)){

  # specify resolution of the trajectory (i.e. granulosity)
  temporal_simplification <- function(df,phi){
  df <- df[df$frame == min(df$frame) | df$frame %% phi  == 0  | df$frame == min(df$frame), ] 
  }

trajectories <- temporal_simplification(trajectories_raw,granulosity[i])

for (j in 1:length(window_size)){

    #make sure that all trajectories have at least as many fixes as window size 
    length <- ddply(trajectories, .(traj), transform, N = length(x))
    trajectories <- trajectories[length$N > window_size[j], ]
    test <- trajectories#[trajectories$traj == 45,]
    
    # rolling diff function
    roll_diff <- function(x) rollapply(x, 2, function(x) diff(x), by.column=F, fill = NA, align = "center")
    
    diff_x <- ave(test$x, test[c("traj")], FUN = roll_diff)
    diff_y <- ave(test$y, test[c("traj")], FUN = roll_diff)
    
    id <- as.data.frame(cbind(test$traj,test$frame, row.names=NULL))
    names(id) <- c("traj","frame")
    disp <- as.data.frame(cbind(id,diff_x,diff_y))
    
    # gross_displacement (works not properly with more than one trajectory)
    disp$disp <- sqrt(disp$diff_x^2+disp$diff_y^2)
    
    # use rollapply for each trajectory
    gd_extract <- function(x) rollapply(x, window_size[j], sum, fill=NA, by.column=F, align = "center")
    gd <- ave(disp$disp, disp[c("traj")], FUN = gd_extract)
    gd <- as.data.frame(cbind(gd,id))
    
    # rolling diff function between first and last observation (to calculate net displacement)
    roll_diff_window <- function(x) rollapply(x, width=window_size[j]+1, function(x) diff(x,(window_size[j]))^2, fill = NA, align="center", by.column=FALSE)
    
    x_net <- ave(test$x, test[c("traj")], FUN = roll_diff_window)
    y_net <- ave(test$y, test[c("traj")], FUN = roll_diff_window)
    
    nd <- sqrt(x_net+y_net)
    nd <- as.data.frame(cbind(nd,id))
    
    if (window_size[j] %% 2 == 0) nd$frame <- nd$frame-granulosity[i]
    
    SI <- merge(gd,nd,by=c("traj","frame"))
    SI <- transform(SI, SI = nd/gd)
    SI <- SI[order(SI$traj,SI$frame), ]
    
    SI$granulosity <- granulosity[i]
    SI$window_size <- window_size[j]

    if (i==1 & j == 1) SI_full <- SI
    if (i>=1 & j>1) SI_full <- rbind(SI_full,SI)
}
}




# function to plot the trajectory and the corresponding MSSI
plot_MSSI <- function(data,random=T,trajectory_select=select_traj){
  
if (random){select_traj <- sample(SI_full$traj,5,replace=F)}

  for (k in 1:length(select_traj)){
  traj <-  subset(SI_full, traj == select_traj[k])

  traj_MSSI <- ggplot(traj, aes(frame, window_size)) +
               geom_tile(data=traj, aes(width=as.numeric(granulosity),height=1,fill = SI)) + 
               scale_fill_gradient2(mid=muted("red"), high="blue", low="white", guide = "colourbar", limits=c(0,1))+
               theme(legend.position="bottom")+
               facet_wrap(~granulosity,ncol=4)

  traj_plot <-  subset(trajectories_raw, traj == select_traj[k])
                gg_traj <- ggplot(traj_plot, aes(x,y,label=frame)) + 
                geom_text(size=3) +
                xlim(0,2048)+
                ylim(0,2048)

  print(grid.arrange(gg_traj,traj_MSSI))

}
}

plot_MSSI(SI_full,random=T)


