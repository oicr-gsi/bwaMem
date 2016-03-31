#!/usr/bin/perl -w

=head2 SYNOPSIS
 
 Using latest file provisioning data this script will check for
 new reports produced with a specified workflow and link to the files
 in appropriate directory

 ./autolinker.pl --study PCSI --data-dir /.mount/labs/PDE/data --token COMP --subdir qc --expand yes --workflow RnaSeqQc

 In case any new reports expanded from a zip archive the script will send an email to the specified
 address with the list of html reports for linking

=cut

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use constant DEBUG=>0;

my $MAIL = "mail";
my($study, $outdir, $datadir, $token, $subdir, $expand, $wf, $latest, $email, @reports);
my $USAGE = "autolinker.pl --study [study id] --data-dir [directory for linking] --token [string to use for filtering] --subdir[subdir to expand to] --workflow [name of workflow] --latest [latest report, gzipped]\n";

my $result = GetOptions("study=s"        => \$study,
                        "data-dir=s"     => \$datadir,
                        "token=s"        => \$token,
                        "subdir=s"       => \$subdir,
                        "expand=s"       => \$expand,
                        "workflow=s"     => \$wf,
                        "email=s"        => \$email,
                        "latest=s"       => \$latest);

if (!$study || !$datadir || !$subdir || !$wf ) { die $USAGE; }
$expand ||="yes";
$datadir =~m!/$! or $datadir.="/";
$latest ||="/.mounts/labs/seqprodbio/private/backups/sqwprod-db.hpc.oicr.on.ca/seqware_files_report_latest.gz";

# =======Checks =======

$latest =~/\.gz$/ or die "Latest report should have .gz extension";
$expand =~/yes|no/ or die "Expand should be set to yes or no";

my $tempfile = "temp_$$";
`zcat $latest | grep $study | grep $token | grep -i $wf | cut -f 2,8,14,47 > $tempfile`;

&validate($tempfile);
&process($tempfile);

# Remove tmp file
`rm $tempfile`;

# email the list of links to reports
if (@reports && @reports > 0) {
 my $message = join("\n",@reports);
 `echo \"$message\" | $MAIL -s "New RNAseq Reports expanded for $study, need to make links on OICR wiki" $email`;
}

=head2

 Validate subroutine makes sure
 we have a well-formed fields in our data file

=cut

sub validate {
 my $f = shift @_;
 open(FILE,"<$f") or die "Couldn't read from file with data";
 while(<FILE>) {
   chomp;
   my @temp = split("\t");
   @temp == 4 or die "Malformed file, won't continue";
   $temp[0] eq $study or die "Line with unrelated study [$temp[0]]  data, aborting";
   $temp[1]=~/^$study/ or die "Malformed donor id [$temp[1]], no study info found";
   $temp[2]=~/^$study/ or die "Unrelated library [$temp[2]] data, aborting";
   
   ($temp[3] && -f $temp[3]) or die "Couldn't validate result file, non-existent file?"
 }
 close FILE;
}

=head2

 Process all entries, create directory if not present
 link to a file

=cut

sub process {
 my $f = shift @_;
 open(FILE,"<$f") or die "Couldn't read from file with data";
 while(<FILE>) {
   chomp;
   my @temp = split("\t");
   my $donor = $temp[1];
   $donor=~s/_//g; #remove all underscores
   my $basedir = $datadir.$donor;
   my $sd = join("/",($basedir,$temp[2],$subdir));
   
   my $filebase = basename($temp[3]);
   if ($filebase =~/fastq/ && $filebase=~/_WG_/) {
     print STDERR "Skipping WG file\n" if DEBUG;
     next;
   } # Hard-coded filter
   if (!-d $sd) {
     print STDERR "Will create directory [$sd]\n" if DEBUG;
     `mkdir -p $sd`;
   } else {
     print STDERR "Directory $sd exists, will not create it\n" if DEBUG;
   }

   my $lf = join("/",($sd,$filebase));
   if (!-e $lf) {
     print STDERR "Will link out to file [$temp[3]] in [$sd]\n" if DEBUG;
     `ln -s $temp[3] -t $sd`;
     if ($lf=~/\.zip$/ && $expand=~/yes/) {
       print STDERR "Will expand [$filebase]\n";
       `cd $sd && unzip $filebase && cd $Bin`;
       my $report = `ls -lh $sd/*/*.html | cut -f 9`;
       if ($report && -e $report) {push (@reports,$report);}
     }
   }
 }
 close FILE;

}
