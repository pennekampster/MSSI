#' Function to calculate the Multiscale Straightness Index
#' 
#' Takes a dataframe with the X and Y coordinates of a trajectory and calculates the MSSI for different window sizes and 
#' granulosities. 
#' @param data A dataframe containing the X and Y coordinates
#' @param uniqueID The unique identifier for each trajectory in the dataframe
#' @param time The column containing the information on the time for each fix
#' @param window_size A numeric value or vector specifying the window sizes over which the MSSI should be calculated
#' @param granulosity A numeric value or vector specifying the temporal resolution of the data before the MSSI is calculated
#' @return A dataframe containing the MSSI for each fix in addition to the input data
#' @import dplyr zoo tcltk
#' @export
#' @references C.M. Postlethwaite, P. Brown, und T.E. Dennis: A New Multi-Scale Measure for Analysing Animal Movement Data. 
#' Journal of Theoretical Biology 317 (2013): 175-185. doi:10.1016/j.jtbi.2012.10.007. 

calculate_MSSI <- function(data,uniqueID="traj",time="frame",window_size,granulosity){

original_id <- eval(parse(text=paste0("data$",uniqueID)))

# sort columns
data <- data[,c(paste0(uniqueID),paste0(time),"X","Y")]

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
trajectories <- as.data.frame(data %>%
                group_by(uniqueID) %>%
                filter(time == min(time) | time %% granulosity[i]  == 0  | time == max(time))%>%
                arrange(uniqueID,time))

for (j in 1:length(window_size)){
  
   Sys.sleep(0.1)
   setTkProgressBar(pb, j, label=paste(round(j/length(window_size)*100, 0), "% done of run ",i, "(out of",length(granulosity),")"))
  
    #make sure that all trajectories have at least as many fixes as window size 
    length <- as.data.frame(trajectories %>%
              group_by(uniqueID) %>%
              mutate(N = length(x)))
   
    trajectories <- trajectories[length$N > window_size[j], ]

    if(nrow(trajectories)==0){ break
    
    } else {
    
    # specify rolling diff function to calculate the displacement between subsequent x or y coordinates
    roll_diff <- function(x) rollapply(x, 2, function(x) diff(x), by.column=F, fill = NA, align = "center")
 
    # run rolling diff function per trajectory
    disp <- as.data.frame(trajectories %>%
                         group_by("uniqueID") %>%
                         transform(diff_x = roll_diff(x), diff_y=roll_diff(y)))
  
    # calculate displacement based on diffs in x and y for subsequent fixes
    disp$disp <- sqrt(disp$diff_x^2+disp$diff_y^2)
    disp <- disp[, c(1,2,7)] 
    
    # use rollapply to sum displacement into gross displacement for each trajectory
    gd_extract <- function(x) rollapply(x, window_size[j], sum, fill=NA, by.column=F, align = "center")

    gd <- as.data.frame(disp %>%
          group_by(uniqueID) %>%
          mutate(gd = gd_extract(disp)))

    gd$disp <- NULL

    # specify rolling diff function between first and last observation (to calculate net displacement)
    roll_diff_window <- function(x) rollapply(x, width=window_size[j]+1, function(x) diff(x,(window_size[j]))^2, fill = NA, align="center", by.column=FALSE)
    
    nd <- as.data.frame(trajectories %>%
          group_by("uniqueID") %>%
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



