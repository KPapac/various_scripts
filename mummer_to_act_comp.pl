#!/usr/bin/perl -w

open (FILE, $ARGV[0]) ||
  die "Can't open $ARGV[0]: $!";

while (<FILE>) {
    chomp;
    s/\s+/ /g;
    s/^\s//;
    s/\s$//;
    @line = split " ", $_;
    if (defined $line[9] && $line[9] =~ /\d+$/) {
	if ($line[9] == "100.00" && $line[1] == $line[4] && $line[0] == $line[3] && $line[11] eq $line[12]){
	}else{
	    print "100 ",$line[9]," ",$line[0]," ",$line[1]," ",$line[11]," ",$line[3]," ",$line[4]," ",$line[12],"\n";
	}
    }
}
