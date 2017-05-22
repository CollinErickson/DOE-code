save_plots_as_gif <- function() {
  # Need to set right directory
  "ffmpeg -f image2 -framerate 1 -i image%d.jpg video.avi"
}