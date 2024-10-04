#!/usr/bin/perl

use strict;
use warnings;

# Get the input file name from the command line argument
my $filename = shift @ARGV;

# Open the input file for reading
open(my $input_fh, '<', $filename) or die "Cannot open $filename: $!";

# Process each line of the input file
while (my $line = <$input_fh>) {
    # Replace tab with " 0 " using the substitute operator
    $line =~ s/\t/ 0 /;

    # Print the modified line
    print $line;
}

# Close the input file
close($input_fh);
