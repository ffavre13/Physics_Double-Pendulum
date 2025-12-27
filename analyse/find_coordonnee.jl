using VideoIO, Images, ImageView

video = VideoIO.openvideo("./analyse/video/First_Video_2s.mp4")
img = read(video)

imshow(img)