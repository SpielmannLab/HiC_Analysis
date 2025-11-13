#!/usr/bin/awk -f
# Script to convert unique.loop from diffpeakachu to .bedpe format for APA
#
BEGIN {
    FS = "\t";
    OFS = "\t";
    print "#chr", "x1", "x2", "chr2" "y1", "y2", "name", "score", "strand1", "strand2", "color";
}

{
    print $0
} 
