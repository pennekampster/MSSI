library(plyr)
library(zoo)
library(ggplot2)
library(gridExtra)

#raw data
trajectories_raw <- read.csv("/Users/Frank/Desktop/video_compress/data/test_videos/3 - trajectory data/ParticleLinker_Data00072.ijout.txt.txt",sep = ",")


# specify window size and granulosity (NB they must be > 3)
granulosity <- seq(1,20,4)
window_size <- 2:20


for (i in 1:length(granulosity)){

#trajectories <- read.csv("/Users/Frank/Dropbox/test_traj.csv",sep = ";")
trajectories <- read.csv("/Users/Frank/Desktop/video_compress/data/test_videos/3 - trajectory data/ParticleLinker_Data00072.ijout.txt.txt",sep = ",")

# specify resolution of the trajectory (i.e. granulosity)
temporal_simplification <- function(df,phi){
  df <- df[df$frame == min(df$frame) | df$frame %% phi  == 0  | df$frame == min(df$frame), ] 
}

trajectories <- temporal_simplification(trajectories,granulosity[i])


for (j in 1:length(window_size)){

#make sure that all trajectories have at least as many fixes as window size 
length <- ddply(trajectories, .(traj), transform, N = length(x))
trajectories <- trajectories[length$N > window_size[j], ]
test <- trajectories#[trajectories$traj == 45,]

# # function working for a single trajectory
# gd <- rollapply(test, width=window_size, function(x) sum(sqrt((diff(x[,3])^2)+(diff(x[,4])^2))), by.column=FALSE, align="left")
# 
# min_frame <- rollapply(test, width=window_size, FUN = function(x){x[which(min(x[,2])==x[,2]),]}, align="left", by.column=FALSE)
# max_frame <- rollapply(test, width=window_size, FUN = function(x){x[which(max(x[,2])==x[,2]),]}, align="left", by.column=FALSE)
# 
# nd <- sqrt((min_frame[,"x"]-max_frame[,"x"])^2+(min_frame[,"y"]-max_frame[,"y"])^2)
# SI <- nd/gd
# 
# SI

# generalization to multiple trajectories
# ave uses factors
#test$traj <- factor(test$traj)
#levels(test$traj)

# rolling diff function
roll_diff <- function(x) rollapply(x, 2, function(x) diff(x), by.column=F, fill = NA, align = "center")

diff_x <- ave(test$x, test[c("traj")], FUN = roll_diff)
diff_y <- ave(test$y, test[c("traj")], FUN = roll_diff)

#diff_x <- na.omit(diff_x)
#diff_y <- na.omit(diff_y)

id <- as.data.frame(cbind(test$traj,test$frame, row.names=NULL))
names(id) <- c("traj","frame")
disp <- as.data.frame(cbind(id,diff_x,diff_y))
#disp <- disp[complete.cases(disp),]

# gross_displacement (works not properly with more than one trajectory)
disp$disp <- sqrt(disp$diff_x^2+disp$diff_y^2)
#gd <- rollapply(disp$disp, window_size-1, sum, by.column=F, align = "left")

# use rollapply for each trajectory
#gd_extract <- function(x) rollapply(x, window_size[j]-1, sum, fill=NA, by.column=F, align = "center")
gd_extract <- function(x) rollapply(x, window_size[j], sum, fill=NA, by.column=F, align = "center")

gd <- ave(disp$disp, disp[c("traj")], FUN = gd_extract)
gd <- as.data.frame(cbind(gd,id))


#plot(test$x,test$y)

# rolling min & max frame functions; returns the minimum frame per trajectory
#get_min_frame <- function(x) rollapply(x, width=window_size, min, fill = NA, align="left", by.column=FALSE)
#get_max_frame <- function(x) rollapply(x, width=window_size, max, fill = NA, align="left", by.column=FALSE)

# rolling diff function between first and last observation (to calculate net displacement)
roll_diff_window <- function(x) rollapply(x, width=window_size[j]+1, function(x) diff(x,(window_size[j]))^2, fill = NA, align="center", by.column=FALSE)

x_net <- ave(test$x, test[c("traj")], FUN = roll_diff_window)
y_net <- ave(test$y, test[c("traj")], FUN = roll_diff_window)

# x_net <- na.omit(x_net) 
# y_net <- na.omit(y_net) 

#x_net <- x_net[!is.na(x_net)]
#y_net <- y_net[!is.na(y_net)]

nd <- sqrt(x_net+y_net)

# nd_data <- as.data.frame(cbind(nd,disp$traj))
# 
# shift <- function(x){
#   x <- x[2:(length(x))]
#   x <- c(x,-999)
#   x[x==-999] <- NA
#   return(x)
# }
# 
# nd <- ave(nd_data$nd, disp[c("traj")], FUN = shift)

nd <- as.data.frame(cbind(nd,id))
#out <- na.omit(nd)/na.omit(gd)

# # shift in nd relative to gd needs to be investigated!!!!!!!!!!!!!!!!
# 
#  if (length(window_size)==1 & window_size[j] == 1) {SI <- nd/gd
#                                                     SI <- as.data.frame(cbind(id,SI), row.names=NULL)}
#                                                     
#                                                     
#  if (length(window_size)>1 & j == 1) {SI <- nd/gd
#                                       SI <- as.data.frame(cbind(id,SI), row.names=NULL)}
#  
#  
#  if (length(window_size)>1 & j > 1) {SI <- nd[2:length(nd)]/gd[1:(length(gd)-1)]
#                                      SI <- as.data.frame(cbind(id[2:length(nd), ],SI), row.names=NULL)}

#if (window_size[j] %% 2 == 0) nd$frame <- nd$frame-1
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




# plot the MSSI
select_traj <- sample(SI_full$traj,20,replace=F)
  
#select_traj <- 1:2
  
for (k in 1:length(select_traj)){

#traj <- subset(SI_full, traj == select_traj[k] & granulosity == 2)
traj <-  subset(SI_full, traj == select_traj[k])
  
  
traj_MSSI <- ggplot(traj, aes(frame, window_size)) +
  geom_tile(data=traj, aes(width=as.numeric(granulosity),height=1,fill = SI)) + 
  scale_fill_gradientn(colours=c("red","violet", "darkblue"),values=c(0,0.5,1),guide="colorbar") +
  theme(legend.position="bottom")+
  facet_wrap(~granulosity)

traj_plot <-  subset(trajectories_raw, traj == select_traj[k])
gg_traj <- ggplot(traj_plot, aes(x,y,label=frame)) + 
geom_text(size=3) +
xlim(0,2048)+
ylim(0,2048)

print(grid.arrange(gg_traj,traj_MSSI))
}




