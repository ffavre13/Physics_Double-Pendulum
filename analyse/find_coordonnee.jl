using VideoIO, Images, ImageView

"""
Open the first frame of the video for finding the pivot coordinate
"""

video = VideoIO.openvideo("./analyse/video/First_Video_2s.mp4")
img = read(video)

imshow(img)