# Calculate a summary image of metaphasic spindles

# The script is adapted to images that are 512 by 512 pixels
# Save the images as jpg gray scale files in ImageJ
# Save control spindle images with the pattern "_ctrl" in the name
# Save the RNase-treated spindle images with the pattern "_RNase" in the name
# Save all images in the same file

# The script analyses first the RNase-treated spindle images and then the control spindle images, based on the name of the images-

# Start RStudio

#### Set working directory ####
# setwd("~/Spindle_Analysis")

# Get packages
library(imager)
library(raster)
library(plotly)

####################
# Print messages 
####################

Prot1 <- readline(prompt="What is the primary protein you are interested in (name): ") # e.g., TUB. Usually the one labelling the microtubules
Prot2 <- readline(prompt="What is the secondary protein you are interested in (name or none): ") # e.g., PLA

# Create file patterns
pattern1 <- as.character(Prot1)
pattern2 <- as.character(Prot2)

if (pattern2 == "none") { # There only one protein to analyse

  ####################
  #### Read files ####
  ####################
  
# Read files from the folder
# Get the number of files with RNase-treated spindles
n_spdl_R <- length(list.files(".", pattern = paste0(pattern1, "_RNase*") ))
print( paste("There are", n_spdl_R, "RNase-treated spindles to analyze", sep = " "))
# Get the number of files with control spindles
n_spdl_C <- length(list.files(".", pattern = paste0(pattern1, "_ctrl*") ))
print( paste("There are", n_spdl_C, "control spindles to analyze", sep = " "))

# Prepare empty variable to store the images
list_im_mean_R <- list()
spindle_sizes_R <- c()

# Print message
print("We start now the analysis of the RNase-treated spindles.")
print("Please, click on the two spindle poles of each spindle.")


###########################################################
#### Click on the poles for the RNase-treated spindles ####
###########################################################

# Make a loop to get the mean spindle size and save the images after shift and rotation
# Get the list of files -> TAKE CARE OF THE PATTERN!
list_files_R <- list.files(".", pattern = paste0(pattern1, "_RNase*") )
for (i in 1:n_spdl_R) {
  
  # Image name
  im_name <- list_files_R[i]
  im_name <- substr(im_name, 1, nchar(im_name)-4)
  
  # Load the images to select the poles
  im_R <- load.image(list_files_R[i])
  plot(im_R)                      # Show the image
  
  # Get the coordinates of the poles
  pole1 <- locator(n = 1)
  pole2 <- locator(n = 1)
  
  # Report the size of the spindle
  size <- sqrt( (pole2$x-pole1$x)^2 + (pole2$y-pole1$y)^2 )
  print( paste("RNase-treated spindle size is:", round(size, digits = 2), "pixels", sep = " "))
  
  # Save the size of the spindles in a vector
  spindle_sizes_R <- c(spindle_sizes_R,size)
  
  # Re-center the image around the spindle
  cntr <- NULL
  cntr$x <- abs(pole2$x - pole1$x)/2 + min(pole1$x, pole2$x) 
  cntr$y <- abs(pole2$y - pole1$y)/2 + min(pole1$y, pole2$y)
  
  # Shift image
  shiftx <- 256 - cntr$x # because the starting image is 512 px by 512 px
  shifty <- 256 - cntr$y
  im_shift_R <- imshift(im_R,shiftx,shifty)
  
  # Determine the axis of rotation to get the spindle aligned vertically       
  rot_axis <- atan( (pole2$x - pole1$x) / (pole2$y - pole1$y) )
  rot_axis <- (rot_axis * 180/pi)   # We don't need to inverse the sign because the y-axis is inverted.
  
  # Rotate the image
  im_rot_R <- imrotate(im_shift_R, rot_axis)
  
  # Save images -> TAKE CARE OF THE NAME! 
  save.image(im_rot_R, paste0(im_name,"_",pattern1,"_RNASE_rotated.jpg"), quality = 1)
}

###########################################
# Repeat the same for the control spindles
###########################################

#### Click on the poles for the control spindles ####
# Empty variables first
#list_im_300_C <- list()
list_im_mean_C <- list()
spindle_sizes_C <- c()
print("We continue now with the analysis of the control spindles.")
print("Please, click on the two spindle poles of each spindle.")
# Make a loop to get the mean spindle size and save the images after shift and rotation
# Get the list of files -> TAKE CARE OF THE PATTERN!
list_files_C <- list.files(".", pattern = paste0(pattern1, "_ctrl*") )
for (i in 1:n_spdl_C) {
  
  # Image name
  im_name <- list_files_C[i]
  im_name <- substr(im_name, 1, nchar(im_name)-4)
  
  # Load the images to select the poles
  im_C <- load.image(list_files_C[i])
  plot(im_C)                      # Show the image
  
  # Get the coordinates of the poles
  pole1 <- locator(n = 1)
  pole2 <- locator(n = 1)
  
  # Report the size of the spindle
  size <- sqrt( (pole2$x-pole1$x)^2 + (pole2$y-pole1$y)^2 )
  print( paste("Control spindle size is:", round(size, digits = 2), "pixels", sep = " "))
  
  # Save the size of the spindles in a vector
  spindle_sizes_C <- c(spindle_sizes_C,size)
  
  # Re-center the image around the spindle
  cntr <- NULL
  cntr$x <- abs(pole2$x - pole1$x)/2 + min(pole1$x, pole2$x) 
  cntr$y <- abs(pole2$y - pole1$y)/2 + min(pole1$y, pole2$y)
  
  # Shift image
  shiftx <- 256 - cntr$x # because the starting image is 512 px by 512 px
  shifty <- 256 - cntr$y
  im_shift_C <- imshift(im_C,shiftx,shifty)
  
  # Determine the axis of rotation to get the spindle aligned vertically       
  rot_axis <- atan( (pole2$x - pole1$x) / (pole2$y - pole1$y) )
  rot_axis <- (rot_axis * 180/pi)   # We don't need to inverse the sign because the y-axis is inverted.
  
  # Rotate the image
  im_rot_C <- imrotate(im_shift_C, rot_axis)
  
  # Save images -> TAKE CARE OF THE NAME!
  save.image(im_rot_C, paste0(im_name,"_",pattern1,"_CTRL_rotated.jpg"), quality = 1)
}

######################################################################################################
# Save the variable which contains the sizes of the individual spindles
save(spindle_sizes_C, file = "spindle_size_C.Rdata")
save(spindle_sizes_R, file = "spindle_size_R.Rdata")
print( paste("Average Control spindle size is:", round(mean(spindle_sizes_C), digits = 0), "pixels", sep = " "))
print( paste("Average RNase spindle size is:", round(mean(spindle_sizes_R), digits = 0), "pixels", sep = " "))


###############################################
# Second loop to calculate the cropped images
###############################################

# RNase-treated spindles
# Get the list of files -> TAKE CARE OF THE PATTERN!
list_files_R <- list.files(".", pattern = paste0(pattern1, "_RNASE_rotated*") )
for (i in 1:n_spdl_R) {
  
  # Load the images to select the poles
  im_rot_R <- load.image(list_files_R[i])
  #plot(im_rot_R)                      # Show the image
  
  # Keep also a version of the spindles that is resized/normalized to the average size of the group
  resize_factor_mean_R <- mean(spindle_sizes_R)/spindle_sizes_R[i] # size in px! Will be the same for one sample (all spindles)
  im_resize_mean_R <- resize(im_rot_R,round(width(im_rot_R)*resize_factor_mean_R),round(height(im_rot_R)*resize_factor_mean_R))
  
  # Crop image with the mean sized spindle in the middle
  cntr_x_mean <- round( width(im_resize_mean_R)/2, digits = 0)
  cntr_y_mean <- round( height(im_resize_mean_R)/2, digits = 0)
  im_crop_mean_R <- imsub(im_resize_mean_R, x %inr% c(cntr_x_mean-150,cntr_x_mean+150), y %inr% c(cntr_y_mean-200,cntr_y_mean+200))

  # Store the image in a list
  #list_im_300_R[[i]] <- im_crop_300_R
  list_im_mean_R[[i]] <- im_crop_mean_R
  
  # Save the images for later
  #save.image(im_crop_300_R, paste("TPX2_R_",i,"_crop_300.jpg", sep =""), quality = 1)
  save.image(im_crop_mean_R, paste0(pattern1,"_R_",i,"_crop_MEAN.jpg"), quality = 1)
  }

# Control spindles
# Get the list of files -> TAKE CARE OF THE PATTERN!
list_files_C <- list.files(".", pattern = "_CTRL_rotated*")
for (i in 1:n_spdl_C) {
  
  # Load the images to select the poles
  im_rot_C <- load.image(list_files_C[i])
  #plot(im_rot_C)                      # Show the image

  # Keep also a version of the spindles that is resized/normalized to the average size of the group
  resize_factor_mean_C <- mean(spindle_sizes_C)/spindle_sizes_C[i] # size in px!
  im_resize_mean_C <- resize(im_rot_C,round(width(im_rot_C)*resize_factor_mean_C),round(height(im_rot_C)*resize_factor_mean_C))
  
  # Crop image with the mean sized spindle in the middle
  cntr_x_mean <- round( width(im_resize_mean_C)/2, digits = 0)
  cntr_y_mean <- round( height(im_resize_mean_C)/2, digits = 0)
  im_crop_mean_C <- imsub(im_resize_mean_C, x %inr% c(cntr_x_mean-150,cntr_x_mean+150), y %inr% c(cntr_y_mean-200,cntr_y_mean+200))
  
  # Store the image in a list
  list_im_mean_C[[i]] <- im_crop_mean_C
  
  # Save the images for later
  save.image(im_crop_mean_C, paste0(pattern1,"_C_",i,"_crop_MEAN.jpg"), quality = 1)
}


#########################################
#### Calculate the average structure ####
#########################################

# Calculate an average image for the TPX2 signal in the RNase-treated sample
im_mat_R <- as.matrix(list_im_mean_R[[1]])
for (i in 2:n_spdl_R) {
  im_mat_R <- im_mat_R + as.matrix( list_im_mean_R[[i]] )
}
im_mat_R <- as.cimg(im_mat_R/n_spdl_R)
save.image(im_mat_R, paste0(pattern1,"_R_","crop_MEAN_spindle.jpg"), quality = 1)

# Calculate an average image for the TPX2 signal in the control sample
im_mat_C <- as.matrix(list_im_mean_C[[1]])
for (i in 2:n_spdl_C) {
  im_mat_C <- im_mat_C + as.matrix( list_im_mean_C[[i]] )
}
im_mat_C <- as.cimg(im_mat_C/n_spdl_C)
save.image(im_mat_C, paste0(pattern1,"_C_","crop_MEAN_spindle.jpg", sep =""), quality = 1)

#### Create montage ####
# Can be done in ImageJ


#############################
#### Plot average images ####
#############################

# Plot average images from the saved file
# CTRL image
plot(load.image( paste0(pattern1,"_C_crop_MEAN_spindle.jpg") ),
     rescale = F,
     axes = F
)
# RNASE image
plot(load.image( paste0(pattern1,"_R_crop_MEAN_spindle.jpg") ),
     rescale = F,
     axes = F
     )


#######################################################
#### Get the intensity profiles for RNase spindles ####
#######################################################

# Get the values of the middle vertical through the average image
mat_R <- matrix(as.vector(im_mat_R), nrow = dim(im_mat_R)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
profile_R <- as.numeric(mat_R[round(dim(im_mat_R)[1]/2),])*255

# Get the values of the middle vertical through all images to add standard error of the mean
# The cropped images are stored in the variables list_im_mean_R
ori_mat_R <- NULL
ori_mat_R <- matrix(as.vector(list_im_mean_R[[1]]), nrow = dim(im_mat_R)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
ori_profile_R <- as.numeric(ori_mat_R[round(dim(im_mat_R)[1]/2),])*255

# Make a dataframe with the profile values ordered per column.
ori_profile_R <- data.frame(ori_profile_R)

# Calculate the standard deviation
for (i in 2:n_spdl_R) {
  ori_mat_R <- matrix(as.vector(list_im_mean_R[[i]]), nrow = dim(im_mat_R)[1], byrow = F)
  ori_profile_R <- cbind( ori_profile_R, data.frame( as.numeric(ori_mat_R[round(dim(im_mat_R)[1]/2),])*255 ) )
}

# Calculate now the standard error of the mean for each point of the profile.
# Add a column to the data frame and then take the values into a vector.
ori_profile_R$SD_M <- apply(ori_profile_R, 1, function(x) {
  round( sd(x), digits = 2)
})

SD_M_R <- as.numeric(ori_profile_R$SD_M)

# Calculate the curves mean +/- SE (SD/sqrt(n_spdl))
SE_p_R <- SD_M_R/sqrt(n_spdl_R) + profile_R
SE_m_R <- profile_R - SD_M_R/sqrt(n_spdl_R)


######################################################
#### Get the intensity profiles for CTRL spindles ####
######################################################

# Get the values of the middle vertical through the average image
mat_C <- matrix(as.vector(im_mat_C), nrow = dim(im_mat_C)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
profile_C <- as.numeric(mat_C[round(dim(im_mat_C)[1]/2),])*255

# Get the values of the middle vertical through all images
# The cropped images are stored in the variables list_im_mean_C
ori_mat_C <- NULL
ori_mat_C <- matrix(as.vector(list_im_mean_C[[1]]), nrow = dim(im_mat_C)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
ori_profile_C <- as.numeric(ori_mat_C[round(dim(im_mat_C)[1]/2),])*255

# Make a dataframe with the profile values ordered per column.
ori_profile_C <- data.frame(ori_profile_C)

# Calculate the standard deviation
for (i in 2:n_spdl_C) {
  ori_mat_C <- matrix(as.vector(list_im_mean_C[[i]]), nrow = dim(im_mat_C)[1], byrow = F)
  ori_profile_C <- cbind( ori_profile_C, data.frame( as.numeric(ori_mat_C[round(dim(im_mat_C)[1]/2),])*255 ) )
}

# Calculate now the standard error of the mean for each point of the profile.
# Add a column to the data frame and then take the values into a vector.
ori_profile_C$SD_M <- apply(ori_profile_C, 1, function(x) {
  round( sd(x), digits = 2)
})

SD_M_C <- as.numeric(ori_profile_C$SD_M)

# Calculate the curves mean +/- SE (SD/sqrt(n_spdl))
SE_p_C <- SD_M_C/sqrt(n_spdl_C) + profile_C
SE_m_C <- profile_C - SD_M_C/sqrt(n_spdl_C)


##############################################################
#### Calculate add the position of the poles on the graph ####
##############################################################

# Average size in pixel for the RNase spindles
RNase_AVG_size <- mean(spindle_sizes_R)
round(RNase_AVG_size)
half_spindle_R <- round(RNase_AVG_size/2)
# Average size in pixel for the CTRL spindles
CTRL_AVG_size <- mean(spindle_sizes_C)
round(CTRL_AVG_size)
half_spindle_C <- round(CTRL_AVG_size/2)
# The poles are positioned equidistant (half-spindle respectively) from the middle of the graph

################################################################################
# save the whole environment here, because all the needed data were produced
base::save.image(file = paste0(pattern1, "_Spindle_Analysis_Environment.Rdata") )
# Can be recall with the function load() when the working directory is properly set.
# load ("Spindle_Analysis_Environment.Rdata")

################################################################################
#### Figure with both profiles on the same plot ####
ylim <- max(profile_C, profile_R) + 10
plot(
  profile_R,
  col = "darkred", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
  lwd = 2,
  type = "l",
  lty = 1,
  xlab = "Position",
  ylab = paste0(pattern1," intensity"),
  ylim = c(0, ylim),
  font.lab = 2,
  main = paste0("Average ",pattern1, " fluorescence distribution")
)
# Add profile for the control spindles
lines(
  profile_C,
  col = "darkgreen", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
  lwd = 2,
  lty = 1)

# Add the standard Error of the mean (SE) for both profiles
x_SE_R <- c(1:length(SE_p_R),length(SE_p_R):1)
y_SE_R <- c(SE_p_R, rev(SE_m_R))
x_SE_C <- c(1:length(SE_p_C),length(SE_p_C):1)
y_SE_C <- c(SE_p_C, rev(SE_m_C))
polygon(x_SE_R, y_SE_R, lty = 6, lwd = 1, col=rgb(0.55, 0.14, 0.14,0.3), border = "darkred")
polygon(x_SE_C, y_SE_C, lty = 6, lwd = 1, col=rgb(0.3, 0.55, 0.45,0.3), border = "darkgreen")

# Add position of the poles as dashed lines
abline(v=c(200-half_spindle_C,200+half_spindle_C), lty = c(2,2), lwd = c(1,1), col = c("darkgreen","darkgreen"))
abline(v=c(200-half_spindle_R,200+half_spindle_R), lty = c(2,2), lwd = c(1,1), col = c("darkred","darkred"))

# Add legend at the top
legend("topright",
       legend = c("Control","RNase-treated"),
       col = c("darkgreen","darkred"), ### HERE SAME COLOR AS ABOVE ##################### 
       lty = 1,
       lwd = 3,
       cex = 0.7)


#### Save graph as image ####
jpeg( paste0("Average_",pattern1,"_Distribution.jpeg") ,width = 8, height = 8, units = 'in', res = 200)

plot(
  profile_R,
  col = "darkred", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
  lwd = 2,
  type = "l",
  lty = 1,
  xlab = "Position",
  ylab = paste0(pattern1," intensity"),
  ylim = c(0,ylim),
  font.lab = 2,
  main = paste0("Average ",pattern1, " fluorescence distribution")
)
# Add profile for the control spindles
lines(
  profile_C,
  col = "darkgreen", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
  lwd = 2,
  lty = 1)

# Add the standard Error of the mean (SE) for both profiles
x_SE_R <- c(1:length(SE_p_R),length(SE_p_R):1)
y_SE_R <- c(SE_p_R, rev(SE_m_R))
x_SE_C <- c(1:length(SE_p_C),length(SE_p_C):1)
y_SE_C <- c(SE_p_C, rev(SE_m_C))
polygon(x_SE_R, y_SE_R, lty = 6, lwd = 1, col=rgb(0.55, 0.14, 0.14,0.3), border = "darkred")
polygon(x_SE_C, y_SE_C, lty = 6, lwd = 1, col=rgb(0.3, 0.55, 0.45,0.3), border = "darkgreen")

# Add position of the poles as dashed lines
abline(v=c(200-half_spindle_C,200+half_spindle_C), lty = c(2,2), lwd = c(1,1), col = c("darkgreen","darkgreen"))
abline(v=c(200-half_spindle_R,200+half_spindle_R), lty = c(2,2), lwd = c(1,1), col = c("darkred","darkred"))

# Add legend at the top
legend("topright",
       legend = c("Control","RNase-treated"),
       col = c("darkgreen","darkred"), ### HERE SAME COLOR AS ABOVE ##################### 
       lty = 1,
       lwd = 3,
       cex = 0.7)

dev.off()


##########################################################################################################
# END OF: if (pattern2 == "none") and beginning of ELSE loop - meaning there is another protein to analyse
##########################################################################################################

} else { 
  
  ####################
  #### Read files ####
  ####################
  
  # Read files from the folder
  # Get the number of files with RNase-treated spindles
  n_spdl_R <- length(list.files(".", pattern = paste0("_RNase_",pattern1) ))
  print( paste("There are", n_spdl_R, "RNase-treated spindles to analyze", sep = " "))
  # Get the number of files with control spindles
  n_spdl_C <- length(list.files(".", pattern = paste0("_ctrl_",pattern1) ))
  print( paste("There are", n_spdl_C, "control spindles to analyze", sep = " "))
  
  # Prepare empty variable to store the images
  list_im_mean_R <- list()
  list_im_mean_R2 <- list()
  spindle_sizes_R <- c()
  
  # Print messages
  print("We start now the analysis of the RNase-treated spindles.")
  print("Please, click on the two spindle poles of each spindle.")
  
  
  ###########################################################
  #### Click on the poles for the RNase-treated spindles ####
  ###########################################################
  
  # Make a loop to get the mean spindle size and save the images after shift and rotation
  # Get the list of files -> TAKE CARE OF THE PATTERN!
  list_files_R1 <- list.files(".", pattern = paste0("_RNase_",pattern1) ) # protein 1
  list_files_R2 <- list.files(".", pattern = paste0("_RNase_",pattern2) ) # protein 2
  for (i in 1:n_spdl_R) {
    
    # Image name
    im_name <- list_files_R1[i]
    im_name <- substr(im_name, 1, nchar(im_name)-4)
    
    # Load the images to select the poles
    im_R <- load.image(list_files_R1[i])
    plot(im_R)                      # Show the image with the protein labeling the microtubules
    
    # Get the coordinates of the poles
    pole1 <- locator(n = 1)
    pole2 <- locator(n = 1)
    
    # Report the size of the spindle
    size <- sqrt( (pole2$x-pole1$x)^2 + (pole2$y-pole1$y)^2 )
    print( paste("RNase-treated spindle size is:", round(size, digits = 2), "pixels", sep = " "))
    
    # Save the size of the spindles in a vector
    spindle_sizes_R <- c(spindle_sizes_R,size)
    
    # Re-center the image around the spindle
    cntr <- NULL
    cntr$x <- abs(pole2$x - pole1$x)/2 + min(pole1$x, pole2$x) 
    cntr$y <- abs(pole2$y - pole1$y)/2 + min(pole1$y, pole2$y)
    
    # Shift image
    shiftx <- 256 - cntr$x # because the starting image is 512 px by 512 px
    shifty <- 256 - cntr$y
    im_shift_R1 <- imshift(im_R,shiftx,shifty)
    
    # Determine the axis of rotation to get the spindle aligned vertically       
    rot_axis <- atan( (pole2$x - pole1$x) / (pole2$y - pole1$y) )
    rot_axis <- (rot_axis * 180/pi)   # We don't need to inverse the sign because the y-axis is inverted.
    
    # Rotate the image
    im_rot_R <- imrotate(im_shift_R1, rot_axis)
    
    # Repeat with the image of protein 2
    im_R2 <- load.image(list_files_R2[i])
    im_name2 <- list_files_R2[i]
    im_name2 <- substr(im_name2, 1, nchar(im_name2)-4)
    im_shift_R2 <- imshift(im_R2,shiftx,shifty)
    im_rot_R2 <- imrotate(im_shift_R2, rot_axis)
    
    # Save images -> TAKE CARE OF THE NAME! 
    save.image(im_rot_R, paste0(im_name,"_",pattern1,"_RNASE_rotated.jpg"), quality = 1)
    save.image(im_rot_R2, paste0(im_name2,"_",pattern2,"_RNASE_rotated.jpg"), quality = 1)
  }
  
  
  ##########################################
  # Repeat the same for the control spindles
  ##########################################
  
  #### Click on the poles for the control spindles ####
  # Empty variables first
  list_im_mean_C <- list()
  list_im_mean_C2 <- list()
  spindle_sizes_C <- c()
  print("We continue now with the analysis of the control spindles.")
  print("Please, click on the two spindle poles of each spindle.")
  # Make a loop to get the mean spindle size and save the images after shift and rotation
  # Get the list of files -> TAKE CARE OF THE PATTERN!
  
  list_files_C1 <- list.files(".", pattern = paste0("_ctrl_",pattern1) ) # protein 1
  list_files_C2 <- list.files(".", pattern = paste0("_ctrl_",pattern2) ) # protein 2
  
  for (i in 1:n_spdl_C) {
    
    # Image name
    im_name <- list_files_C1[i]
    im_name <- substr(im_name, 1, nchar(im_name)-4)
    
    # Load the images to select the poles
    im_C <- load.image(list_files_C1[i])
    plot(im_C)                      # Show the image
    
    # Get the coordinates of the poles
    pole1 <- locator(n = 1)
    pole2 <- locator(n = 1)
    
    # Report the size of the spindle
    size <- sqrt( (pole2$x-pole1$x)^2 + (pole2$y-pole1$y)^2 )
    print( paste("Control spindle size is:", round(size, digits = 2), "pixels", sep = " "))
    
    # Save the size of the spindles in a vector
    spindle_sizes_C <- c(spindle_sizes_C,size)
    
    # Re-center the image around the spindle
    cntr <- NULL
    cntr$x <- abs(pole2$x - pole1$x)/2 + min(pole1$x, pole2$x) 
    cntr$y <- abs(pole2$y - pole1$y)/2 + min(pole1$y, pole2$y)
    
    # Shift image
    shiftx <- 256 - cntr$x # because the starting image is 512 px by 512 px
    shifty <- 256 - cntr$y
    im_shift_C <- imshift(im_C,shiftx,shifty)
    
    # Determine the axis of rotation to get the spindle aligned vertically       
    rot_axis <- atan( (pole2$x - pole1$x) / (pole2$y - pole1$y) )
    rot_axis <- (rot_axis * 180/pi)   # We don't need to inverse the sign because the y-axis is inverted.
    
    # Rotate the image
    im_rot_C <- imrotate(im_shift_C, rot_axis)
    
    # Repeat with the image of protein 2
    im_C2 <- load.image(list_files_C2[i])
    im_name2 <- list_files_C2[i]
    im_name2 <- substr(im_name2, 1, nchar(im_name2)-4)
    im_shift_C2 <- imshift(im_C2,shiftx,shifty)
    im_rot_C2 <- imrotate(im_shift_C2, rot_axis)
    
    # Save images -> TAKE CARE OF THE NAME!
    save.image(im_rot_C, paste0(im_name,"_",pattern1,"_CTRL_rotated.jpg"), quality = 1)
    save.image(im_rot_C2, paste0(im_name2,"_",pattern2,"_CTRL_rotated.jpg"), quality = 1)
  }
  
  
  ######################################################################################################
  # Save the variable which contains the sizes of the individual spindles
  save(spindle_sizes_C, file = "spindle_size_C.Rdata")
  save(spindle_sizes_R, file = "spindle_size_R.Rdata")
  print( paste("Average Control spindle size is:", round(mean(spindle_sizes_C), digits = 0), "pixels", sep = " "))
  print( paste("Average RNase spindle size is:", round(mean(spindle_sizes_R), digits = 0), "pixels", sep = " "))
  
  
  #############################################
  # Second loop to calculate the cropped images
  #############################################
  
  # RNase-treated spindles
  # Get the list of files -> TAKE CARE OF THE PATTERN!
  list_files_R1 <- list.files(".", pattern = paste0(pattern1, "_RNASE_rotated*") )
  list_files_R2 <- list.files(".", pattern = paste0(pattern2, "_RNASE_rotated*") )
  for (i in 1:n_spdl_R) {
    
    # Load the images to select the poles
    im_rot_R <- load.image(list_files_R1[i])
    im_rot_R2 <- load.image(list_files_R2[i])
    #plot(im_rot_R)                      # Show the image

    
    # Keep also a version of the spindles that is resized/normalized to the average size of the group
    resize_factor_mean_R <- mean(spindle_sizes_R)/spindle_sizes_R[i] # size in px! Will be the same for one sample (all spindles)
    im_resize_mean_R <- resize(im_rot_R,round(width(im_rot_R)*resize_factor_mean_R),round(height(im_rot_R)*resize_factor_mean_R))
    im_resize_mean_R2 <- resize(im_rot_R2,round(width(im_rot_R2)*resize_factor_mean_R),round(height(im_rot_R2)*resize_factor_mean_R))
    
    # Crop image with the mean sized spindle in the middle
    cntr_x_mean <- round( width(im_resize_mean_R)/2, digits = 0)
    cntr_y_mean <- round( height(im_resize_mean_R)/2, digits = 0)
    im_crop_mean_R <- imsub(im_resize_mean_R, x %inr% c(cntr_x_mean-150,cntr_x_mean+150), y %inr% c(cntr_y_mean-200,cntr_y_mean+200))
    im_crop_mean_R2 <- imsub(im_resize_mean_R2, x %inr% c(cntr_x_mean-150,cntr_x_mean+150), y %inr% c(cntr_y_mean-200,cntr_y_mean+200))
    
    # Store the image in a list
    list_im_mean_R[[i]] <- im_crop_mean_R
    list_im_mean_R2[[i]] <- im_crop_mean_R2
    
    # Save the images for later
    save.image(im_crop_mean_R, paste0(pattern1,"_R_",i,"_crop_MEAN.jpg"), quality = 1)
    save.image(im_crop_mean_R2, paste0(pattern2,"_R_",i,"_crop_MEAN.jpg"), quality = 1)    
  }
  
  # Control spindles
  # Get the list of files -> TAKE CARE OF THE PATTERN!
  list_files_C1 <- list.files(".", pattern = paste0(pattern1, "_CTRL_rotated*") )
  list_files_C2 <- list.files(".", pattern = paste0(pattern2, "_CTRL_rotated*") )
  
  for (i in 1:n_spdl_C) {
    
    # Load the images to select the poles
    im_rot_C <- load.image(list_files_C1[i])
    im_rot_C2 <- load.image(list_files_C2[i])
    #plot(im_rot_C)                      # Show the image
    
    # Keep also a version of the spindles that is resized/normalized to the average size of the group
    resize_factor_mean_C <- mean(spindle_sizes_C)/spindle_sizes_C[i] # size in px!
    im_resize_mean_C <- resize(im_rot_C,round(width(im_rot_C)*resize_factor_mean_C),round(height(im_rot_C)*resize_factor_mean_C))
    im_resize_mean_C2 <- resize(im_rot_C2,round(width(im_rot_C2)*resize_factor_mean_C),round(height(im_rot_C2)*resize_factor_mean_C))
    
    # Crop image with the mean sized spindle in the middle
    cntr_x_mean <- round( width(im_resize_mean_C)/2, digits = 0)
    cntr_y_mean <- round( height(im_resize_mean_C)/2, digits = 0)
    im_crop_mean_C <- imsub(im_resize_mean_C, x %inr% c(cntr_x_mean-150,cntr_x_mean+150), y %inr% c(cntr_y_mean-200,cntr_y_mean+200))
    im_crop_mean_C2 <- imsub(im_resize_mean_C2, x %inr% c(cntr_x_mean-150,cntr_x_mean+150), y %inr% c(cntr_y_mean-200,cntr_y_mean+200))
    
    # Store the image in a list
    list_im_mean_C[[i]] <- im_crop_mean_C
    list_im_mean_C2[[i]] <- im_crop_mean_C2
    
    # Save the images for later
    save.image(im_crop_mean_C, paste0(pattern1,"_C_",i,"_crop_MEAN.jpg"), quality = 1)
    save.image(im_crop_mean_C2, paste0(pattern2,"_C_",i,"_crop_MEAN.jpg"), quality = 1)
  }
  
  
  #########################################
  #### Calculate the average structure ####
  #########################################
  
  # Calculate an average image for the protein signal in the RNase-treated sample
  im_mat_R <- as.matrix(list_im_mean_R[[1]])
  im_mat_R2 <- as.matrix(list_im_mean_R2[[1]])
  for (i in 2:n_spdl_R) {
    im_mat_R <- im_mat_R + as.matrix( list_im_mean_R[[i]] )
    im_mat_R2 <- im_mat_R2 + as.matrix( list_im_mean_R2[[i]] )
  }
  im_mat_R <- as.cimg(im_mat_R/n_spdl_R)
  save.image(im_mat_R, paste0(pattern1,"_R_","crop_MEAN_spindle.jpg"), quality = 1)
  im_mat_R2 <- as.cimg(im_mat_R2/n_spdl_R)
  save.image(im_mat_R2, paste0(pattern2,"_R_","crop_MEAN_spindle.jpg"), quality = 1)
  
  # Calculate an average image for the protein signal in the control sample
  im_mat_C <- as.matrix(list_im_mean_C[[1]])
  im_mat_C2 <- as.matrix(list_im_mean_C2[[1]])
  for (i in 2:n_spdl_C) {
    im_mat_C <- im_mat_C + as.matrix( list_im_mean_C[[i]] )
    im_mat_C2 <- im_mat_C2 + as.matrix( list_im_mean_C2[[i]] )
  }
  im_mat_C <- as.cimg(im_mat_C/n_spdl_C)
  save.image(im_mat_C, paste0(pattern1,"_C_","crop_MEAN_spindle.jpg"), quality = 1)
  im_mat_C2 <- as.cimg(im_mat_C2/n_spdl_C)
  save.image(im_mat_C2, paste0(pattern2,"_C_","crop_MEAN_spindle.jpg"), quality = 1)
  
  
  
  #########################################
  #### Calculate the standard deviation  ##
  #########################################
  
  # list of images
  # list_im_mean_C, list_im_mean_C1, list_im_mean_R, list_im_mean_R2
  
  # Standard deviation of the lists
  SD_C <- parsd(list_im_mean_C)
  SD_C2 <- parsd(list_im_mean_C2)
  SD_R <- parsd(list_im_mean_R)
  SD_R2 <- parsd(list_im_mean_R2)
  
  plot(SD_C)
  plot(SD_C2)
  plot(SD_R)
  plot(SD_R2)

  save.image(  SD_C, paste0(pattern1,"_C_","crop_SD_spindle.jpg"), quality = 1)
  save.image(  SD_C2, paste0(pattern1,"_C2_","crop_SD_spindle.jpg"), quality = 1)
  save.image(  SD_R, paste0(pattern1,"_R_","crop_SD_spindle.jpg"), quality = 1)  
  save.image(  SD_R2, paste0(pattern1,"_R2_","crop_SD_spindle.jpg"), quality = 1)
  
  # Same with the variance
  plot(parvar(list_im_mean_C2))
  plot(parvar(list_im_mean_R2))  
  
  save.image(  parsd(list_im_mean_C2), paste0(pattern1,"_VarC2_","crop_SD_spindle.jpg"), quality = 1)
  save.image(  parvar(list_im_mean_R2), paste0(pattern1,"_VarR2_","crop_SD_spindle.jpg"), quality = 1)
  
  
  p_im_C2 <- plot_ly(z = ~as.matrix(parsd(list_im_mean_C2)),
                  colorscale = list(c(0,1),
                                    c("rgb(100,100,100)","rgb(139,0,0)")
                  )) %>% 
    add_surface() %>% layout(title = "AVG PLA signal")
  p_im_C2
  
  p_im_R2 <- plot_ly(z = ~as.matrix(parsd(list_im_mean_R2)),
                     colorscale = list(c(0,1),
                                       c("rgb(100,100,100)","rgb(139,0,0)")
                     )) %>% 
    add_surface() %>% layout(title = "AVG PLA signal")
  p_im_R2
  
  #### Create montage ####
  # Can be done in ImageJ
  
  #############################
  #### Plot average images ####
  #############################
  
  # Plot average images from the saved file
  # CTRL mean images
  plot(load.image( paste0(pattern1,"_C_crop_MEAN_spindle.jpg") ),
       rescale = F,
       axes = F
  )
  plot(load.image( paste0(pattern2,"_C_crop_MEAN_spindle.jpg") ),
       rescale = F,
       axes = F
  )
  # RNASE mean images
  plot(load.image( paste0(pattern1,"_R_crop_MEAN_spindle.jpg") ),
       rescale = F,
       axes = F
  )
  plot(load.image( paste0(pattern2,"_R_crop_MEAN_spindle.jpg") ),
       rescale = F,
       axes = F
  )

  
  #######################################################
  #### Get the intensity profiles for RNase spindles ####
  #######################################################

  # Get the values of the middle vertical through the average image
  mat_R <- matrix(as.vector(im_mat_R), nrow = dim(im_mat_R)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  profile_R <- as.numeric(mat_R[round(dim(im_mat_R)[1]/2),])*255
  mat_R2 <- matrix(as.vector(im_mat_R2), nrow = dim(im_mat_R2)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  profile_R2 <- as.numeric(mat_R2[round(dim(im_mat_R2)[1]/2),])*255
  
  # Get the values of the middle vertical through all images to add standard error of the mean
  # The cropped images are stored in the variables list_im_mean_R
  ori_mat_R <- NULL
  ori_mat_R <- matrix(as.vector(list_im_mean_R[[1]]), nrow = dim(im_mat_R)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  ori_profile_R <- as.numeric(ori_mat_R[round(dim(im_mat_R)[1]/2),])*255
  ori_mat_R2 <- NULL
  ori_mat_R2 <- matrix(as.vector(list_im_mean_R2[[1]]), nrow = dim(im_mat_R2)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  ori_profile_R2 <- as.numeric(ori_mat_R2[round(dim(im_mat_R2)[1]/2),])*255
  
  # Make a dataframe with the profile values ordered per column.
  ori_profile_R <- data.frame(ori_profile_R)
  ori_profile_R2 <- data.frame(ori_profile_R2)
  
  # Calculate the standard deviation
  for (i in 2:n_spdl_R) {
    ori_mat_R <- matrix(as.vector(list_im_mean_R[[i]]), nrow = dim(im_mat_R)[1], byrow = F)
    ori_profile_R <- cbind( ori_profile_R, data.frame( as.numeric(ori_mat_R[round(dim(im_mat_R)[1]/2),])*255 ) )
  }
  for (i in 2:n_spdl_R) {
    ori_mat_R2 <- matrix(as.vector(list_im_mean_R2[[i]]), nrow = dim(im_mat_R2)[1], byrow = F)
    ori_profile_R2 <- cbind( ori_profile_R2, data.frame( as.numeric(ori_mat_R2[round(dim(im_mat_R2)[1]/2),])*255 ) )
  }
  
  # Calculate now the standard error of the mean for each point of the profile.
  # Add a column to the data frame and then take the values into a vector.
  ori_profile_R$SD_M <- apply(ori_profile_R, 1, function(x) {
    round( sd(x), digits = 2)
  })
  ori_profile_R2$SD_M <- apply(ori_profile_R2, 1, function(x) {
    round( sd(x), digits = 2)
  })
  
  SD_M_R <- as.numeric(ori_profile_R$SD_M)
  SD_M_R2 <- as.numeric(ori_profile_R2$SD_M)
  
  # Calculate the curves mean +/- SE (SD/sqrt(n_spdl))
  SE_p_R <- SD_M_R/sqrt(n_spdl_R) + profile_R
  SE_m_R <- profile_R - SD_M_R/sqrt(n_spdl_R)
  SE_p_R2 <- SD_M_R2/sqrt(n_spdl_R) + profile_R2
  SE_m_R2 <- profile_R2 - SD_M_R2/sqrt(n_spdl_R)
  
  
  #######################################################
  #### Get the intensity profiles for CTRL spindles ####
  #######################################################
  
  # Get the values of the middle vertical through the average image
  mat_C <- matrix(as.vector(im_mat_C), nrow = dim(im_mat_C)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  profile_C <- as.numeric(mat_C[round(dim(im_mat_C)[1]/2),])*255
  mat_C2 <- matrix(as.vector(im_mat_C2), nrow = dim(im_mat_C2)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  profile_C2 <- as.numeric(mat_C2[round(dim(im_mat_C2)[1]/2),])*255
  
  # Get the values of the middle vertical through all images
  # The cropped images are stored in the variables list_im_mean_C
  ori_mat_C <- NULL
  ori_mat_C <- matrix(as.vector(list_im_mean_C[[1]]), nrow = dim(im_mat_C)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  ori_profile_C <- as.numeric(ori_mat_C[round(dim(im_mat_C)[1]/2),])*255
  ori_mat_C2 <- NULL
  ori_mat_C2 <- matrix(as.vector(list_im_mean_C2[[1]]), nrow = dim(im_mat_C2)[1], byrow = F) # !!! The matrix is rotated -90 degrees.
  ori_profile_C2 <- as.numeric(ori_mat_C2[round(dim(im_mat_C2)[1]/2),])*255
  
  # Make a dataframe with the profile values ordered per column.
  ori_profile_C <- data.frame(ori_profile_C)
  ori_profile_C2 <- data.frame(ori_profile_C2)
  
  # Calculate the standard deviation
  for (i in 2:n_spdl_C) {
    ori_mat_C <- matrix(as.vector(list_im_mean_C[[i]]), nrow = dim(im_mat_C)[1], byrow = F)
    ori_profile_C <- cbind( ori_profile_C, data.frame( as.numeric(ori_mat_C[round(dim(im_mat_C)[1]/2),])*255 ) )
  }
  for (i in 2:n_spdl_C) {
    ori_mat_C2 <- matrix(as.vector(list_im_mean_C2[[i]]), nrow = dim(im_mat_C2)[1], byrow = F)
    ori_profile_C2 <- cbind( ori_profile_C2, data.frame( as.numeric(ori_mat_C2[round(dim(im_mat_C2)[1]/2),])*255 ) )
  }
  
  # Calculate now the standard error of the mean for each point of the profile.
  # Add a column to the data frame and then take the values into a vector.
  ori_profile_C$SD_M <- apply(ori_profile_C, 1, function(x) {
    round( sd(x), digits = 2)
  })
  ori_profile_C2$SD_M <- apply(ori_profile_C2, 1, function(x) {
    round( sd(x), digits = 2)
  })
  
  SD_M_C <- as.numeric(ori_profile_C$SD_M)
  SD_M_C2 <- as.numeric(ori_profile_C2$SD_M)
  
  # Calculate the curves mean +/- SE (SD/sqrt(n_spdl))
  SE_p_C <- SD_M_C/sqrt(n_spdl_C) + profile_C
  SE_m_C <- profile_C - SD_M_C/sqrt(n_spdl_C)
  SE_p_C2 <- SD_M_C2/sqrt(n_spdl_C) + profile_C2
  SE_m_C2 <- profile_C2 - SD_M_C2/sqrt(n_spdl_C)
  
  
  ##############################################################
  #### Calculate add the position of the poles on the graph ####
  ##############################################################
  
  # Average size in pixel for the RNase spindles
  RNase_AVG_size <- mean(spindle_sizes_R)
  round(RNase_AVG_size)
  half_spindle_R <- round(RNase_AVG_size/2)
  # Average size in pixel for the CTRL spindles
  CTRL_AVG_size <- mean(spindle_sizes_C)
  round(CTRL_AVG_size)
  half_spindle_C <- round(CTRL_AVG_size/2)
  # The poles are positioned equidistant (half-spindle respectively) from the middle of the graph
  
  
  ################################################################################
  # save the whole environment here, because all the needed data were produced
  ################################################################################
  
  base::save.image(file = paste0(pattern1,"_",pattern2, "_Spindle_Analysis_Environment.Rdata") )
  # Can be recall with the function load() when the working directory is properly set.
  # load ("XXX_Analysis_Environment.Rdata")
  
  
  ################################################################################
  #### Figure for protein 1 with both profiles on the same plot ####
  ################################################################################
  
  ylim <- max(profile_C, profile_R) + 10
  plot(
    profile_R,
    col = "darkred", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    type = "l",
    lty = 1,
    xlab = "Position",
    ylab = paste0(pattern1," intensity"),
    ylim = c(0, ylim),
    font.lab = 2,
    main = paste0("Average ",pattern1, " fluorescence distribution")
  )
  # Add profile for the control spindles
  lines(
    profile_C,
    col = "darkgreen", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    lty = 1)
  
  # Add the standard Error of the mean (SE) for both profiles
  x_SE_R <- c(1:length(SE_p_R),length(SE_p_R):1)
  y_SE_R <- c(SE_p_R, rev(SE_m_R))
  x_SE_C <- c(1:length(SE_p_C),length(SE_p_C):1)
  y_SE_C <- c(SE_p_C, rev(SE_m_C))
  polygon(x_SE_R, y_SE_R, lty = 6, lwd = 1, col=rgb(0.55, 0.14, 0.14,0.3), border = "darkred")
  polygon(x_SE_C, y_SE_C, lty = 6, lwd = 1, col=rgb(0.3, 0.55, 0.45,0.3), border = "darkgreen")
  
  # Add position of the poles as dashed lines
  abline(v=c(200-half_spindle_C,200+half_spindle_C), lty = c(2,2), lwd = c(1,1), col = c("darkgreen","darkgreen"))
  abline(v=c(200-half_spindle_R,200+half_spindle_R), lty = c(2,2), lwd = c(1,1), col = c("darkred","darkred"))
  
  # Add legend at the top
  legend("topright",
         legend = c("Control","RNase-treated"),
         col = c("darkgreen","darkred"), ### HERE SAME COLOR AS ABOVE ##################### 
         lty = 1,
         lwd = 3,
         cex = 0.7)
  
  
  #### Save graph as image ####
  jpeg( paste0("Average_",pattern1,"_Distribution.jpeg") ,width = 8, height = 8, units = 'in', res = 200)
  
  plot(
    profile_R,
    col = "darkred", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    type = "l",
    lty = 1,
    xlab = "Position",
    ylab = paste0(pattern1," intensity"),
    ylim = c(0,ylim),
    font.lab = 2,
    main = paste0("Average ",pattern1, " fluorescence distribution")
  )
  # Add profile for the control spindles
  lines(
    profile_C,
    col = "darkgreen", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    lty = 1)
  
  # Add the standard Error of the mean (SE) for both profiles
  x_SE_R <- c(1:length(SE_p_R),length(SE_p_R):1)
  y_SE_R <- c(SE_p_R, rev(SE_m_R))
  x_SE_C <- c(1:length(SE_p_C),length(SE_p_C):1)
  y_SE_C <- c(SE_p_C, rev(SE_m_C))
  polygon(x_SE_R, y_SE_R, lty = 6, lwd = 1, col=rgb(0.55, 0.14, 0.14,0.3), border = "darkred")
  polygon(x_SE_C, y_SE_C, lty = 6, lwd = 1, col=rgb(0.3, 0.55, 0.45,0.3), border = "darkgreen")
  
  # Add position of the poles as dashed lines
  abline(v=c(200-half_spindle_C,200+half_spindle_C), lty = c(2,2), lwd = c(1,1), col = c("darkgreen","darkgreen"))
  abline(v=c(200-half_spindle_R,200+half_spindle_R), lty = c(2,2), lwd = c(1,1), col = c("darkred","darkred"))
  
  # Add legend at the top
  legend("topright",
         legend = c("Control","RNase-treated"),
         col = c("darkgreen","darkred"), ### HERE SAME COLOR AS ABOVE ##################### 
         lty = 1,
         lwd = 3,
         cex = 0.7)
  
  dev.off()
  
  
  ################################################################################
  #### Figure for protein 2 with both profiles on the same plot ####
  ################################################################################
  
  ylim <- max(profile_C2, profile_R2) + 10
  plot(
    profile_R2,
    col = "darkred", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    type = "l",
    lty = 1,
    xlab = "Position",
    ylab = paste0(pattern2," intensity"),
    ylim = c(0, ylim),
    font.lab = 2,
    main = paste0("Average ",pattern2, " fluorescence distribution")
  )
  # Add profile for the control spindles
  lines(
    profile_C2,
    col = "darkgreen", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    lty = 1)
  
  # Add the standard Error of the mean (SE) for both profiles
  x_SE_R <- c(1:length(SE_p_R2),length(SE_p_R2):1)
  y_SE_R <- c(SE_p_R2, rev(SE_m_R2))
  x_SE_C <- c(1:length(SE_p_C2),length(SE_p_C2):1)
  y_SE_C <- c(SE_p_C2, rev(SE_m_C2))
  polygon(x_SE_R, y_SE_R, lty = 6, lwd = 1, col=rgb(0.55, 0.14, 0.14,0.3), border = "darkred")
  polygon(x_SE_C, y_SE_C, lty = 6, lwd = 1, col=rgb(0.3, 0.55, 0.45,0.3), border = "darkgreen")
  
  # Add position of the poles as dashed lines
  abline(v=c(200-half_spindle_C,200+half_spindle_C), lty = c(2,2), lwd = c(1,1), col = c("darkgreen","darkgreen"))
  abline(v=c(200-half_spindle_R,200+half_spindle_R), lty = c(2,2), lwd = c(1,1), col = c("darkred","darkred"))
  
  # Add legend at the top
  legend("topright",
         legend = c("Control","RNase-treated"),
         col = c("darkgreen","darkred"), ### HERE SAME COLOR AS ABOVE ##################### 
         lty = 1,
         lwd = 3,
         cex = 0.7)
  
  
  #### Save graph as image ####
  jpeg( paste0("Average_",pattern2,"_Distribution.jpeg") ,width = 8, height = 8, units = 'in', res = 200)
  
  plot(
    profile_R2,
    col = "darkred", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    type = "l",
    lty = 1,
    xlab = "Position",
    ylab = paste0(pattern2," intensity"),
    ylim = c(0,ylim),
    font.lab = 2,
    main = paste0("Average ",pattern2, " fluorescence distribution")
  )
  # Add profile for the control spindles
  lines(
    profile_C2,
    col = "darkgreen", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
    lwd = 2,
    lty = 1)
  
  # Add the standard Error of the mean (SE) for both profiles
  x_SE_R <- c(1:length(SE_p_R2),length(SE_p_R2):1)
  y_SE_R <- c(SE_p_R2, rev(SE_m_R2))
  x_SE_C <- c(1:length(SE_p_C2),length(SE_p_C2):1)
  y_SE_C <- c(SE_p_C2, rev(SE_m_C2))
  polygon(x_SE_R, y_SE_R, lty = 6, lwd = 1, col=rgb(0.55, 0.14, 0.14,0.3), border = "darkred")
  polygon(x_SE_C, y_SE_C, lty = 6, lwd = 1, col=rgb(0.3, 0.55, 0.45,0.3), border = "darkgreen")
  
  # Add position of the poles as dashed lines
  abline(v=c(200-half_spindle_C,200+half_spindle_C), lty = c(2,2), lwd = c(1,1), col = c("darkgreen","darkgreen"))
  abline(v=c(200-half_spindle_R,200+half_spindle_R), lty = c(2,2), lwd = c(1,1), col = c("darkred","darkred"))
  
  # Add legend at the top
  legend("topright",
         legend = c("Control","RNase-treated"),
         col = c("darkgreen","darkred"), ### HERE SAME COLOR AS ABOVE ##################### 
         lty = 1,
         lwd = 3,
         cex = 0.7)
  
  dev.off()
  
} # END of ELSE loop



############################################################################################
# Individual graphs if needed
############################################################################################

#### Plot for the RNase sample ####
#plot(
#  profile_R,
#  col = "darkred", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
#  lwd = 2,
#  type = "l",
#  xlab = "Position",
#  ylab = paste0(pattern1," intensity"),
#  ylim = c(0,max(profile_R)+10),
#  font.lab = 2,
#  main = paste0("Average ",pattern1, " fluorescence distribution")
#)
#legend("topright",
#       legend = c("RNase"), # TAKE CARE OF THE NAME OF THE PROTEIN FOR THE LEGEND #####################
#       col = c("darkred"), ### HERE SAME COLOR AS ABOVE ##################### 
#       lty = 1,
#       lwd = 2,
#       cex = 0.6,
#       xpd = T)

#Plot the SE on the graph
# use the function polygone
#x_SE_R <- c(1:length(SE_p_R),length(SE_p_R):1)
#y_SE_R <- c(SE_p_R, rev(SE_m_R))
#polygon(x_SE_R, y_SE_R, lty = 6, lwd = 1, border = "black")
# polygon(x_SE, y_SE, col = rgb(250, 142, 130, max = 255), lty = 3, lwd = 1, border = "darkred")

# Add position of the poles as dashed lines
#abline(v=c(200-half_spindle_R,200+half_spindle_R), lty = c(2,2), lwd = c(1,1), col = c("darkred","darkred"))

#############################################
# To draw transparent rectangles on the graph
#rect_transp <- rgb(190, 190, 190, alpha = 130, maxColorValue=255)
#rect(48,0,52,1, density = 100, col = rect_transp, border = NA)
#rect(148,0,152,1, density = 100, col = rect_transp, border = NA)

#### Plot for the control sample ####
#plot(
#  profile_C,
#  col = "darkgreen", ### HERE YOU CAN USE ANY COLOR YOU LIKE #####################
#  lwd = 2,
#  type = "l",
#  xlab = "Position",
#  ylab = paste0(pattern1," intensity"),
#  ylim = c(0,max(profile_C)+10),
#  font.lab = 2,
#  main = paste0("Average ",pattern1, " fluorescence distribution")
#)
#legend("topright",
#       legend = c("CTRL"), # TAKE CARE OF THE NAME OF THE PROTEIN FOR THE LEGEND #####################
#       col = c("darkgreen"), ### HERE SAME COLOR AS ABOVE ##################### 
#       lty = 1,
#       lwd = 2,
#       cex = 0.6,
#       xpd = T)

#Plot the SE on the graph
# use the function polygone
#x_SE_C <- c(1:length(SE_p_C),length(SE_p_C):1)
#y_SE_C <- c(SE_p_C, rev(SE_m_C))
#polygon(x_SE_C, y_SE_C, lty = 6, lwd = 1, border = "black")

# Add position of the poles as dashed lines
#abline(v=c(200-half_spindle_C,200+half_spindle_C), lty = c(2,2), lwd = c(1,1), col = c("darkgreen","darkgreen"))

### Save profiles as table ###
write.csv(ori_profile_C,"ori_profile_WT.csv",col.names = TRUE)
write.csv(ori_profile_R,"ori_profile_KO.csv",col.names = TRUE)