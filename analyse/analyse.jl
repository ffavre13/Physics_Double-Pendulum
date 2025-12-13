using VideoIO

io = VideoIO.open("First_Video_2s.mp4")
f = VideoIO.openvideo(io)
img = read(f)

while !eof(f)
    read!(f, img)
end
close(f)