ifile="$1.md"
ofile="$1.html"
title=$2
pandoc --mathjax -s --metadata title="$title" $ifile -o $ofile
