ffmpeg -framerate 10 -i images/frame_%03d.png -c:v libx264 -pix_fmt yuv420p out.mp4