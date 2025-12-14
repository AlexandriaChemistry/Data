#!/bin/sh

plotxvg -f minus.xvg -yframe 9 -legend_x 0.69 -legend_y 0.99 -xmax 3 -sfx -0.1 -sfy 0.96 -save minus.pdf -noshow -lfs 48 -tickfs 36 -alfs 48 -mksize 16

plotxvg -f plus.xvg -yframe 9 -legend_x 0.69 -legend_y 0.99 -xmax 3 -sfx -0.1 -sfy 0.96  -save plus.pdf -noshow -lfs 48 -tickfs 36 -alfs 48 -mksize 16
