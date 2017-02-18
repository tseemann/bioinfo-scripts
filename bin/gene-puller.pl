#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Bio::SearchIO;
use File::Spec;
use Fatal;

my(@Options, $verbose, $gfn, $dir, $evalue, $tophit);
setOptions();


$dir or err("Please supply a results folder with --dir");
-d $dir and err("Folder '$dir' already exists, sorry.");
msg("Making working folder: $dir");
mkdir $dir or err("Unable to make folder '$dir', sorry.");

$gfn or err("Please supply a input gene FASTA file --genes");
-r $gfn or err("Can't read '$gfn' genes file.");
msg("Gene file: $gfn");

my %cfn;
my $LCP = longest_common_prefix(@ARGV);
for my $cfn (@ARGV) {
  if (-r $cfn) {
#    my $id = id_from_path($cfn);
    my $id = $cfn;
    $id =~ s/$LCP//;
    $id =~ s{/}{_}g;
    msg("Contig file #".(scalar keys %cfn)." aka '$id': $cfn");
    $cfn{$cfn} = $id;
  }
  else {
    err("File '$cfn' is not readable");    
  }  
}
scalar(keys %cfn) or err("Not enough valid contig files were supplied.");
#print STDERR Dumper(\%cfn);

for my $exe (qw(blastall fastacmd formatdb date muscle)) {
  require_exe($exe) ? msg("Found '$exe' - ok") : err("Can't find '$exe' in \$PATH");
}


msg("Creating BLAST databases");
for my $cfn (keys %cfn) {
  my $id = $cfn{$cfn};
  my $cmd = "zcat -f '$cfn' | formatdb -t '$id' -i stdin -n '$dir/$id' -p F -l '$dir/formatdb.log' -V T -o T";
  run_cmd($cmd);
}

msg("BLASTing genes");
my $gio = Bio::SeqIO->new(-file=>$gfn, -format=>'fasta');
while (my $gene = $gio->next_seq) {
  my $tool = $gene->alphabet eq 'dna' ? 'blastn' : 'tblastn';
  my $id = id_from_path( $gene->display_id );
  my $qout = Bio::SeqIO->new(-file=>">$dir/$id.qry", -format=>'fasta');
  $qout->write_seq($gene);
  my $hout = Bio::SeqIO->new(-file=>">$dir/$id.fasta", -format=>'fasta');
  $hout->write_seq($gene);
DATABASE:
  for my $db (values %cfn) {
    msg("Finding $id in $db using $tool ...");
    my $vb = $tophit ? 1 : 99;
#    my $WS = $tool eq 'blastn' ? 5 : 3;
    my $cmd = "blastall -p $tool -i '$dir/$id.qry' -d '$dir/$db' -e $evalue -F F -v $vb -b $vb";
    msg("Running: $cmd");
    open BLS, '-|', $cmd;
    my $bls = Bio::SearchIO->new(-fh=>\*BLS, -format=>'blast');
    my $count=0;
    while (my $res = $bls->next_result) {
      while (my $hit = $res->next_hit) {
        while (my $hsp = $hit->next_hsp) {
	  $count++;
          msg("Got hit #$count to ", $hit->name, " ", $hsp->length);
          my $hs = uc $hsp->hit_string;
          $hs =~ s/-//g;
          my $seq = Bio::Seq->new(
            -id => $db.($count > 1 ? ":hit$count" : ''),
	    -desc => $hit->name.' '.$hsp->start('hit').' '.$hsp->end('hit')
	             .' '.($hsp->strand('hit') > 0 ? '+' : '-'),
	    -seq => $hs,
	  );
          $hout->write_seq($seq);
	  next DATABASE if $tophit;
	}
      }
    }
  }
  # Align them
  msg("Aligning all the $id sequences");
  my $cmd = "muscle -in '$dir/$id.fasta' -out '$dir/$id.aln' -clwstrict -quiet -loga '$dir/muscle.log'";
  run_cmd($cmd);
}

#print STDERR Dumper(\%cfn);
for my $db (values %cfn) {
  msg("deleting temporary BLAST databases for $db");
  for my $ext (qw(nsq nhr nin nsd nsi)) {
    unlink "$dir/$db.$ext";
  }
}
msg("Done! Check the '$dir' folder for your results.");


#-------------------------------------------------------------------------

sub longest_common_prefix {
  my $prefix = shift;
  for (@_) {
    chop $prefix while (! /^\Q$prefix\E/);
  }
  return $prefix;
}

#-------------------------------------------------------------------------

sub id_from_path {
  my($s) = @_;
  my @s = split m{/}, $s;
  $s = $s[-1];
  return $s;
}

#-------------------------------------------------------------------------

sub num_cpu {
  if ($^O =~ m/linux/i) {
    my($num) = qx(grep -c ^processor /proc/cpuinfo);
    chomp $num;
    return $num if $num =~ m/^\d+/;
  }
  return 1;
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
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"g|genes=s",  VAR=>\$gfn, DEFAULT=>'', DESC=>"File name of query gene(s) - ok to mix DNA and PROT"},
    {OPT=>"d|dir=s",  VAR=>\$dir, DEFAULT=>'', DESC=>"Folder to put results in"},
    {OPT=>"e|evalue=f",  VAR=>\$evalue, DEFAULT=>1E-6, DESC=>"e-value cutoff for BLAST similarity"},
    {OPT=>"t|tophit!",  VAR=>\$tophit, DEFAULT=>0, DESC=>"Only include top hit per sample"},
  );

  #(!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print STDERR "Synopsis:\n  Finds, collates and aligns orthologous genes in sets of contigs.\n";
  print STDERR "Usage:\n  $0 [options] -g <genes.fa> -d <outdir> <contigs1.fna [contigs2.fna ...]>\n";
  print STDERR "Options:\n";
  foreach (@Options) {
    printf STDERR "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
