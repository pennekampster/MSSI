#' Function to visualize the MSSI along a trajectory 
#' 
#' Takes the dataframe returned by the calculate_MSSI function and the raw trajectory data to construct the plot
#' @param raw_traj A dataframe containing the X and Y coordinates of the trajectories on which the MSSI was calculated
#' @param data The dataframe returned by the calculate MSSI_function containing the MSSI index for each fix and different window sizes 
#' and granulosities
#' @param uniqueID Column name of the unique identifier for each trajectory in the dataframe
#' @param time Column name containing the information on the time for each fix
#' @param granulosity_choosen A numeric value of granulosity for plotting corresponding to the granulosity values choosen
#' @param random Logical; whether the trajectories to be plotted are randomly sampled from the data 
#' @param N_traj A numeric value specifying the number of trajectories to be plotted
#' @return A plot showing the MSSI for each fix for different window sizes
#' @param trajectory_select Vector containing the unique IDs of the trajectories to be plotted
#' @import ggplot2 gridExtra 
#' @export

plot_MSSI <- function (raw_traj, data, uniqueID = "traj", time = "frame", 
                       granulosity_choosen = 1, random = T, N_traj = 10, trajectory_select = select_traj) 
{
  colnames(raw_traj)[tolower(colnames(raw_traj)) == "x"] <- "x"
  colnames(raw_traj)[tolower(colnames(raw_traj)) == "y"] <- "y"
  colnames(data)[colnames(data) == paste(uniqueID)] <- "uniqueID"
  colnames(data)[colnames(data) == paste(time)] <- "time"
  colnames(raw_traj)[colnames(raw_traj) == paste(uniqueID)] <- "uniqueID"
  colnames(raw_traj)[colnames(raw_traj) == paste(time)] <- "time"
  if (random) {
    select_traj <- sample(data$uniqueID, N_traj, replace = F)
  }
  for (k in 1:length(select_traj)) {
      traj <- subset(data, uniqueID == select_traj[k])
      traj_MSSI <- ggplot(traj, aes(time, window_size)) + 
                 geom_tile(data = subset(traj, granulosity == granulosity_choosen), aes(width = as.numeric(granulosity), height = as.numeric(2), fill = SI)) + 
      scale_fill_gradientn(colours = c("red", "cyan", "black"), guide = "colourbar", limits = c(0, 1)) + 
      ylab("Window size") + 
      xlab("Time") + 
      theme(legend.position = "bottom")
    
      traj_plot <- subset(raw_traj, uniqueID == select_traj[k])
      
      gg_traj <- ggplot(traj_plot, aes(x, y, label = time)) + 
      geom_text(size = 3) +
      ggtitle(paste0("Multiscale straightness index for trajectory: ", select_traj[k]))      
      
      grid.arrange(gg_traj, traj_MSSI, ncol=1, nrow=2)
    
  }
}



