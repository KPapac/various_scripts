#! /usr/bin/perl -w

open (TAB, $ARGV[0]) or die "Can't open $ARGV[0]: $!";
open (FASTA, $ARGV[1]) or die "Can't open $ARGV[1]: $!";

#$found=0;
$fasta=$ARGV[1];
 
while (<FASTA>) {
  chomp;
  if (/^\>/) {
  }else{
      s/\s//g;
      tr/ATGC/atgc/;
      $genome .= $_;
  }
}

$length = length ($genome);

$/="\n";
$i=0;
while (<TAB>) {
  chomp;
  if ($i == 0){
      print "LOCUS       selected             ",$length," bp    dna     circular   UNK\n";
      print "DEFINITION  ",$fasta,"\n";
      print "ACCESSION   unknown\n";
      print "FEATURES             Location/Qualifiers\n"; 
  }
  $i++;
  if (/^FT(.+)/){
      print "  ",$1,"\n";
  }
}
print "ORIGIN\n";
@sequence_chunks = ($genome =~ /.{1,10}/g);
$count=1;
$count2=1;
for $chunk (@sequence_chunks) {
    $count2=1 if $count2 == 7;
    if ($count<10){
	print "        ",$count," ",$chunk," " if $count2 == 1; 
	print $chunk," " if $count2 < 6 && $count2 > 1; 
	print $chunk,"\n" if $count2 == 6;
    }elsif ($count<100){
	print "       ",$count," ",$chunk," " if $count2 == 1; 
	print $chunk," " if $count2 < 6 && $count2 > 1; 
	print $chunk,"\n" if $count2 == 6;
    }elsif ($count<1000){
	print "      ",$count," ",$chunk," " if $count2 == 1; 
	print $chunk," " if $count2 < 6 && $count2 > 1; 
	print $chunk,"\n" if $count2 == 6;
    }elsif ($count<10000){
	print "     ",$count," ",$chunk," " if $count2 == 1; 
	print $chunk," " if $count2 < 6 && $count2 > 1; 
	print $chunk,"\n" if $count2 == 6;
    }elsif ($count<10000){
	print "    ",$count," ",$chunk," " if $count2 == 1; 
	print $chunk," " if $count2 < 6 && $count2 > 1; 
	print $chunk,"\n" if $count2 == 6;
    }elsif ($count<100000){
	print "   ",$count," ",$chunk," " if $count2 == 1; 
	print $chunk," " if $count2 < 6 && $count2 > 1; 
	print $chunk,"\n" if $count2 == 6;
    }elsif ($count<1000000){
	print "  ",$count," ",$chunk," " if $count2 == 1; 
	print $chunk," " if $count2 < 6 && $count2 > 1; 
	print $chunk,"\n" if $count2 == 6;
    }
    $count+=10;
    $count2++ if $count2 < 7;
} 
print "\n";
