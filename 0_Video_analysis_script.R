######################################################################
# R script for analysing video files with BEMOVI (www.bemovi.info)
rm(list=ls())
# load package
library(devtools)
install_github("efronhofer/bemovi", ref="experimental")
library(bemovi)

######################################################################
# VIDEO PARAMETERS

# video frame rate (in frames per second)
fps <- 25
# length of video (in frames)
total_frames <- 500

# measured volume (in microliter)
measured_volume <- 34.4 # for Leica M205 C with 1.6 fold magnification, 
#Sample  height 0.5 mm and Hamamatsu Orca Flash 4

# size of a pixel (in micrometer)
pixel_to_scale <- 4.05 # for Leica M205 C with 1.6 fold magnification, 
#Sample  height 0.5 mm and Hamamatsu Orca Flash 4

# specify video file format (one of "avi","cxd","mov","tiff")
# bemovi only works with avi and cxd. other formats are reformated 
# to avi below
video.format <- "cxd"

# setup
difference.lag <- 10
thresholds <- c(10,255) # don't change the second value
#thresholds <- c(50,255)

######################################################################
# FILTERING PARAMETERS 
# min and max size: area in pixels
particle_min_size <- 5
particle_max_size <- 1000

# number of adjacent frames to be considered for linking particles
trajectory_link_range <- 3
# maximum distance a particle can move between two frames
trajectory_displacement <- 16

# these values are in the units defined by the parameters above: 
# fps (seconds),
#measured_volume (microliters) and pixel_to_scale (micometers)
filter_min_net_disp <- 25
filter_min_duration <- 1
filter_detection_freq <- 0.1
filter_median_step_length <- 3

######################################################################
# MORE PARAMETERS (USUALLY NOT CHANGED)

# set paths to ImageJ and particle linker standalone
IJ.path <- "/home/felix/bin/ImageJ"
to.particlelinker <- "/home/felix/bin/ParticleLinker"

# directories and file names
to.data <- paste(getwd(),"/",sep="")
video.description.folder <- "0_video_description/"
video.description.file <- "video_description.txt"
raw.video.folder <- "1_raw/"
particle.data.folder <- "2_particle_data/"
trajectory.data.folder <- "3_trajectory_data/"
temp.overlay.folder <- "4a_temp_overlays/"
overlay.folder <- "4_overlays/"
merged.data.folder <- "5_merged_data/"
ijmacs.folder <- "ijmacs/"

# RAM allocation
memory.alloc <- c(60000)

# RAM per particle linker instance
memory.alloc.perLinker <- c(10000)

######################################################################
# VIDEO ANALYSIS

# identify particles
locate_and_measure_particles(to.data, raw.video.folder, 
                             particle.data.folder,
                             difference.lag, thresholds, min_size = particle_min_size, 
                             max_size = particle_max_size, IJ.path, memory.alloc)

# link the particles
link_particles(to.data, particle.data.folder, trajectory.data.folder, 
               linkrange = trajectory_link_range, disp = trajectory_displacement, 
               start_vid = 1, memory = memory.alloc, 
               memory_per_linkerProcess = memory.alloc.perLinker)

# merge info from description file and data
merge_data(to.data, particle.data.folder, trajectory.data.folder,
           video.description.folder, video.description.file, merged.data.folder)

# load the merged data
load(paste0(to.data, merged.data.folder, "Master.RData"))

# filter data: minimum net displacement, their duration, the detection 
#frequency and the median step length
trajectory.data.filtered <- filter_data(trajectory.data, 
                                        filter_min_net_disp, filter_min_duration, filter_detection_freq, 
                                        filter_median_step_length)

# summarize trajectory data to individual-based data
morph_mvt <- summarize_trajectories(trajectory.data.filtered, 
                                    calculate.median=F, write = T, to.data, merged.data.folder)

# get Sample  level info
summarize_populations(trajectory.data.filtered, morph_mvt, 
                      write=T, to.data, merged.data.folder, video.description.folder, 
                      video.description.file, total_frames)

# create overlays for validation
create_overlays(trajectory.data.filtered, to.data, 
                merged.data.folder, raw.video.folder, temp.overlay.folder, 
                overlay.folder, 2048, 2048, difference.lag, type = "label", 
                predict_spec = F, IJ.path, contrast.enhancement = 1, 
                memory = memory.alloc)