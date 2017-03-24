#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Spec;
use Fatal;

my(@Options, $verbose, $prot);
setOptions();

require_exe("backtranseq");

my $out = Bio::SeqIO->new(-fh=>\*STDOUT, -format=>'fasta');
for my $fasta (@ARGV) {
  msg("Back translating: $fasta");
  open my $back, '-|', "backtranseq -auto -filter -verbose < $fasta";
  my $in = Bio::SeqIO->new(-fh=>$back, -format=>'fasta');
  while (my $seq = $in->next_seq) {
    msg("Creating 3 test cases for", $seq->id);
    # negative
    my $bg = random_dna( 3 * $seq->length + int(rand(3)) );
    $out->write_seq( Bio::Seq->new(-id=>$seq->id, -desc=>'FALSE', -seq=>$bg) );
    # normal
    substr $bg, $seq->length, $seq->length, $seq->seq;
    $out->write_seq( Bio::Seq->new(-id=>$seq->id, -desc=>'TRUE', -seq=>$bg) );
    # revcom
    substr $bg, $seq->length, $seq->length, $seq->revcom->seq;
    $out->write_seq( Bio::Seq->new(-id=>$seq->id, -desc=>'TRUE', -seq=>$bg) );
  }
}

#-------------------------------------------------------------------------

sub random_dna {
  my($L) = @_;
  my @src = ('A','T','G','C');
  $L or err("random_dna: need length");
  return join( '', map { $src[int(rand(4))] } (1 .. $L) ); 
}

#-------------------------------------------------------------------------

sub run_cmd {
  msg("Running: @_");
  if (system(@_)) {
    msg("ERROR $? : $!");
    exit $?;
  }
}

#-------------------------------------------------------------------------

sub require_exe {
  my($bin) = shift;
  for my $dir (File::Spec->path) {
    my $exe = File::Spec->catfile($dir, $bin);
    return $exe if -x $exe; 
  }
  return;
}
#----------------------------------------------------------------------

sub msg {
  my($time) = qx(date +"%F %T");
  chomp $time;
  print STDERR "[$time] @_\n";
}

#----------------------------------------------------------------------

sub err {
  print STDERR "ERROR: @_\n";
  exit -1;
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
#    {OPT=>"prot!",  VAR=>\$prot, DEFAULT=>0, DESC=>"Make source in protein space"},
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
  print STDERR "Synopsis:\n  Make a test set from target FASTA files\n";
  print STDERR "Usage:\n  $0 [options] [--prot] <genes.faa | genes.fna>\n";
  print STDERR "Options:\n";
  foreach (@Options) {
    printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
