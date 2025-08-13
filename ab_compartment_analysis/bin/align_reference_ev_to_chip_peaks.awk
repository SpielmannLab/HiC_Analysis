#!/usr/bin/env -S awk -f

# Now calculate whether to flip the peaks or not
BEGIN {
  FS = "\t";
  OFS = "\t"
} 
pass == 1 {
  if ($4 == "A") {
    ++n_abins[$1]
    npeaks_in_a[$1] += $7
  } else {
    ++n_bbins[$1]
    npeaks_in_b[$1] += $7
  }
}
pass == 2 {
  if ( !($1 in n_abins) || !($1 in n_bbins) || !(n_bbins[$1] ~ /[[:digit:]]+/) || !(n_bbins[$1] ~ /[[:digit:]]+/) || ((npeaks_in_a[$1] + npeaks_in_b[$1]) == 0) ) {
    printf "The chrom %s does not have either A or B bin or peaks in either of them to align the eigen vector. Skipping this chrom\n", $1 > "/dev/stderr"
    # Choosing not to print unaligned vector elements
    # avg_no_peaks_a[$1] = 0
    # avg_no_peaks_b[$1] = 0
    # $5 = 0
    # $4 = "-"
    # print $0
    next
  }
  avg_no_peaks_a[$1] = npeaks_in_a[$1]/n_abins[$1]
  avg_no_peaks_b[$1] = npeaks_in_b[$1]/n_bbins[$1]
  if ( avg_no_peaks_a[$1] < avg_no_peaks_b[$1] ) {
    $5 = -$5
    if ($4 == "B") { $4 = "A"} else { $4 = "B" }
  }
  print $0
}
END {
  format = "%-10s %-10s %-13s %-16s %-10s %-13s %-16s %-10s\n"
  printf format, "chrom", "n_abins", "npeaks_in_a", "avg_no_peaks_a", "n_bbins", "npeaks_in_b", "avg_no_peaks_b", "inverted?"  > "/dev/stderr"
  for (key in n_abins) {
    if (key in n_bbins) {
      printf format, key, n_abins[key], npeaks_in_a[key], avg_no_peaks_a[key], n_bbins[key], npeaks_in_b[key], avg_no_peaks_b[key], (avg_no_peaks_a[key] < avg_no_peaks_b[key]) > "/dev/stderr"
    n_abins_genome =+ n_abins[key]
    n_peaks_in_a_bins_genome =+ npeaks_in_a[key]
    n_bbins_genome =+ n_bbins[key]
    n_peaks_in_b_bins_genome =+ npeaks_in_b[key]
    }
  } 

  print "" > "/dev/stderr"
  print "GENOME LEVEL STATS" > "/dev/stderr"
  print "------ A compartment ------" > "/dev/stderr"
  printf "Total number of bins in genome: %d\n", n_abins_genome > "/dev/stderr"
  printf "Total number of peaks in genome: %d\n", n_peaks_in_a_bins_genome > "/dev/stderr"
  printf "Average number of peaks per bin in genome: %d\n", n_peaks_in_a_bins_genome/n_abins_genome > "/dev/stderr"
  print "" > "/dev/stderr"
  print "------ B compartment ------" > "/dev/stderr"
  printf "Total number of bins in genome: %d\n", n_bbins_genome > "/dev/stderr"
  printf "Total number of peaks in genome: %d\n", n_peaks_in_b_bins_genome > "/dev/stderr"
  printf "Average number of peaks per bin in genome: %d\n", n_peaks_in_b_bins_genome/n_bbins_genome > "/dev/stderr"
}
