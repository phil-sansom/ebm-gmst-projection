######################################
## Figure 2 - Energy Balance Models ##
######################################

width  = (210/25.4 - 2)/2
height = width*sqrt(2)
pdf("./fig/figure2.pdf", width = width, height = height, pointsize = 12)
par(ann = FALSE, mar = c(1,1,1,1), xaxs = "i", yaxs = "i")

plot (1, 1, type = "n", xlim = c(0,9), ylim = c(0,10), axes = FALSE)

arrows(4, 10, 4, 9, length = 0.1)
text(3.5, 9.5, "N", adj = 1)

segments(0, 9, 8, 9)

arrows(2, 9, 2, 8, length = 0.1)
text(1.5, 8.5, "F", adj = 1)

arrows(6, 8, 6, 9, length = 0.1)
text(5.5, 8.5, expression(paste(plain(k)[1],plain(T)[1])), adj = 1)

rect(0, 7, 8, 8, col = "lightgrey")
text(4, 7.5, expression(paste(plain(C)[1],plain(T)[1])))

arrows(4, 7, 4, 6, length = 0.1)
text(3.5, 6.5, expression(paste(plain(k)[2],"(",plain(T)[1]-plain(T)[2],")")), 
     adj = 1)

rect(0, 4, 8, 6, col = "lightgrey")
text(4, 5, expression(paste(plain(C)[2],plain(T)[2])))

arrows(4, 4, 4, 3, length = 0.1)
text(3.5, 3.5, expression(paste(epsilon,plain(k)[3],
                                "(",plain(T)[2]-plain(T)[3],")")), adj = 1)

rect(0, 0, 8, 3, col = "lightgrey")
text(4, 1.5, expression(paste(plain(C)[3],plain(T)[3])))

segments(7, 9, 7, 8.5)
segments(7, 8.5, 9, 8.5)
segments(9, 8.5, 9, 1.5)
arrows(9, 1.5, 8, 1.5, length = 0.1)
text(8.5, 5, expression(paste("(",1-epsilon,")",plain(k)[3],
                              "(",plain(T)[2]-plain(T)[3],")")), srt = 90)

dev.off()
