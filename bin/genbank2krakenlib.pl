#!/usr/bin/env perl
use strict;

my(@Options, $debug);
setOptions();

for my $infile (@ARGV) {
  print STDERR "Opening: $infile\n";
  my $gi = '';
  my $acc = '';
  my $org = '';
  my $def = '';
  my $in_seq = 0;
  my $dna = '';
  open my $infh, '-|', "gzip -c -d \Q$infile\E";
  while (<$infh>) {
    chomp;
    if (m{^//}) {
      print ">gi|$gi|gb|$acc| $def [$org]\n$dna";
      $in_seq = 0;
      $dna = '';
      next;
    }
    elsif (m/^ORIGIN/) {
      $in_seq = 1;
      next;
    }
    
    if ($in_seq) {
      my $s = substr $_, 10;
      $s =~ s/\s//g;
      $dna .= uc($s);
      $dna .= "\n";
    }
    else {
      if (m/^VERSION.*?GI:(\d+)/) {
        $gi = $1;
      }
      elsif (m/^SOURCE\s+(.*)$/) {
        $org = $1;
      }
      elsif (m/^LOCUS\s+(\S+)/) {
        $acc = $1;
      }
      elsif (m/^DEFINITION\s+(.*)$/) {
        $def = $1;
      }
    }
  }
}
print STDERR "Done\n";

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] < file.gbk > file.fna\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
