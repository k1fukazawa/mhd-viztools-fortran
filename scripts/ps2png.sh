#!/bin/bash
# convert PostScript from slice_plot / polar_map to PNG / PDF (needs ghostscript)
# landscape pages (slice_plot) are rotated automatically; portrait (polar_map) is kept
#   ./ps2png.sh slice.ps [resolution dpi=150]   -> slice-001.png, slice-002.png, ...
#   ./ps2png.sh slice.ps pdf               -> slice.pdf
set -e
ps=$1; res=${2:-150}; base=${ps%.ps}
# landscape if " 90.0 ro" appears in the page setup at the beginning
if head -c 20000 "$ps" | grep -q " 90.0 ro"; then
  ori=(-c "<</Orientation 3>> setpagedevice" -f)
else
  ori=()
fi
if [ "$res" = "pdf" ]; then
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sPAPERSIZE=a4 -dAutoRotatePages=/None \
     -sOutputFile="$base.pdf" "${ori[@]}" "$ps"
  echo "$base.pdf"
else
  gs -q -dNOPAUSE -dBATCH -sDEVICE=png16m -sPAPERSIZE=a4 -r"$res" \
     -dTextAlphaBits=4 -dGraphicsAlphaBits=4 \
     -sOutputFile="$base-%03d.png" "${ori[@]}" "$ps"
  ls "$base"-*.png
fi
