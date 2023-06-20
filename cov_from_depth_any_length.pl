#! /usr/bin/perl -w

open (FILE, $ARGV[0]) or
  die "Can't open $ARGV[0]: $!";

$window_size=$ARGV[1];

while (<FILE>) {
  chomp;
  @line=split /\t/,$_;
  push @pos,$line[1];
  $cov{$line[1]}=$line[2];
}

if (defined $ARGV[2]) {
    $length = $ARGV[2];
}else{
    $length =(pop @pos)+2*$window_size;
} 

$k=0;
for ($i=$window_size;$i<=$length+$window_size;$i+=$window_size){
    $tot_cov = 0;
    foreach $location (@pos) {
	shift @pos if $location < $k;
	last if $location > $i;
	next if $location < $k;
	$tot_cov+=$cov{$location} if $location <= $i && $location >= $k;
    }
    $average_cov=$tot_cov/$window_size;
    print $i,"\t",$average_cov,"\n";
    $i=$i;
    $k=$i;
}
close FILE;
