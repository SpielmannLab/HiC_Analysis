#!/usr/bin/env -S awk -f
# Perform three passes. 
# First pass to check query_ev alignment to ref_ev
# Second pass to flip the query_ev as needed
# Third pass to flip the domain assignments in the query as needed
BEGIN {
  FS = "\t";
  OFS = "\t"
} 
pass == 1 {
  sum_mag[$1] += ($5+$7)^2;
  diff_mag[$1] += ($5-$7)^2;
} 
# The ref compartment that has not been corrected because of lack of chip-signal has been assigned "-". For these, assign the query also to be "-".
(pass == 2) && ($4 != "-") {
  if (diff_mag[$1] > sum_mag[$1]) {
    $7 = -$7 
    if ($6 == "B") { $6 = "A"} else { $6 = "B" }
    }
  print $1,$2,$3,$6,$7,"." >> "query_ev_corrected.bed"
}
# Chosing to not print uncorrected vector elements
# (pass == 2) && ($4 == "-") {
#     $7 = 0
#     $6 = "-"
#   print $1,$2,$3,"-",0,"." >> "query_ev_corrected.bed"
# }
(pass == 3) && ($1 in sum_mag) {
  if (diff_mag[$1] > sum_mag[$1]) {
    $5 = -$5 
    if ($4 == "B") { $4 = "A"} else { $4 = "B" }
    }
  print $0 >> "query_domains_corrected.bed"
}
END {
  for (key in sum_mag) {
    if (diff_mag[key] > sum_mag[key]) {
      printf "Inverted eigen vector of chrom %s\n", key > "correction_stats.txt"
    } 
  } 
}
