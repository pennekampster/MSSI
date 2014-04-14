\name{calculate_MSSI}
\alias{calculate_MSSI}
\title{Function to calculate the Multiscale Straightness Index}
\usage{
calculate_MSSI(data, uniqueID = "traj", time = "frame", window_size,
  granulosity)
}
\arguments{
  \item{data}{A dataframe containing the X and Y
  coordinates}

  \item{uniqueID}{The unique identifier for each trajectory
  in the dataframe}

  \item{time}{The column containing the information on the
  time for each fix}

  \item{window_size}{A numeric value or vector specifying
  the window sizes over which the MSSI should be
  calculated}

  \item{granulosity}{A numeric value or vector specifying
  the temporal resolution of the data before the MSSI is
  calculated}
}
\value{
A dataframe containing the MSSI for each fix in addition to
the input data
}
\description{
Takes a dataframe with the X and Y coordinates of a
trajectory and calculates the MSSI for different window
sizes and granulosities
}

