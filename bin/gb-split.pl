#!/usr/bin/env perl
use strict;

my(@Options, $debug, $template);
setOptions();

my $count = 0;
my $record = '';
my $acc = '';
my $locus = '';

while (<ARGV>) {
  if (m{^//}) {
    $count++;
    $acc ||= "UNKNOWN_$count";
    my $fname = $template;
    $fname =~ s/%i/$count/g;
    $fname =~ s/%a/$acc/g;
    print STDERR "[$count] <$fname> $locus";
    open my $out, '>', $fname;
    print $out $record, '//', "\n";
    $record = '';
    $acc = '';
    $locus = '';
  }
  elsif (m/^LOCUS\s+(\S+)/) {
    $acc = $1;
    $locus = $_;
  }
  $record .= $_;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"debug!",  VAR=>\$debug, DEFAULT=>0, DESC=>"Debug info"},
    {OPT=>"template=s",  VAR=>\$template, DEFAULT=>'%a.gbk', DESC=>"Output template: %a=accession %i=counter"},
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
  print "Usage: $0 [options] file.gbk\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
