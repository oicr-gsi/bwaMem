#!/usr/bin/perl
use warnings;
use strict;

use JSON;
use Data::Dumper;
use Getopt::Std;
use jsonReporter;
use jsonReporter::ConfigData;
use constant DEBUG=>0;
use constant RUNNING=>3;
use constant JSN=>"text/json";
use constant XML=>"text/xml";

$| = 1;

# OOZIE WEBSERVICE:
my @Status = ("FAILED","KILLED","PREP","RUNNING","SUSPENDED","SUCCEEDED","OK");
my $webservice  = jsonReporter::ConfigData->config('webservice');
my $web_basedir = "11000/oozie/v1/job/";
my $web_bulkdir = "11000/oozie/v1/jobs?len=";
my $web_info    = "show=info";
my $web_defi    = "show=definition";

my $recentCutoff = 7 * 24 * 60 * 60; # week
my $recentCutoffTime = "week";
# query the seqware metadb
my $username = jsonReporter::ConfigData->config('username');
my $password = jsonReporter::ConfigData->config('password');
my $dbhost   = jsonReporter::ConfigData->config('dbhost');
my $dbname   = jsonReporter::ConfigData->config('dbname');
my $devmode = 0;
my $timemode= 0;
my $outfile;
my $outdir;


# Control development output with -d option
my %opts = ();
getopts('sdt', \%opts);
if (defined $opts{d}) {
  $devmode = 1;
}

if (defined $opts{t}) {
  $timemode = 1;
}

if (exists $ARGV[0])
{
  if ($ARGV[0] eq "month")
  {
    $recentCutoffTime = $ARGV[0];
    $recentCutoff = 31 * 24 * 60 * 60;
    $outfile = "month.json";
  }
  elsif ($ARGV[0] eq "year")
  {
    $recentCutoffTime = $ARGV[0];
    $recentCutoff = 365 * 24 * 60 * 60;
    $outfile = "year.json";
  }
  elsif ($ARGV[0] eq "decade")
  {
    $recentCutoffTime = $ARGV[0];
    $recentCutoff = 10 * 365 * 24 * 60 * 60;
    $outfile = "decade.json";
  }
  elsif ($ARGV[0] ne "week")
  {
    die "Please run as jsonReport.pl [week|month|year|decade] [output_dir]\n";
  }
  # defaults are fine
}

if (exists $ARGV[1] && -d $ARGV[1]) {
  $outdir = $ARGV[1];
} else {
  die "Need an output directory for writing into";
}

$outfile ||= "week.json";
warn "Connecting to meta_db\n";


my $jr = new jsonReporter($dbhost,$dbname,$username,$password,$webservice,$devmode,$timemode);
my %results = %{$jr->getSWData()};

# Need to know progress
if (defined $results{pending} && scalar(@{$results{pending}}) > 0) {
 warn "Have pending workflows, will parse ".scalar(@{$results{pending}})." runs";
 RUN:
 foreach my $run (@{$results{pending}}) {
    if (!defined $run->{status} || $run->{status} ne $Status[RUNNING]) {
     next RUN;
    }
          # FAKE STATUS IF NO PROGRESS AVAILABLE
          my $progress = "".int rand(100); 
          if (defined $run->{wrun_id} && $run->{wrun_id}=~/oozie/) {
             # Create a request
             my $url = join(":",($webservice,$web_basedir)).$run->{wrun_id}."?";
             my $steps = $jr->countActions($url.$web_defi,XML);
             my $done  = $jr->countActions($url.$web_info,JSN);
             $progress = int($done/$steps*100);
             print STDERR "Got real Progress value $progress\n" if DEBUG; 
          }
          $run->{progress} = $progress;
  } 
}

my $json = JSON->new();
$json = encode_json \%results;

# Atomic update
$outdir.="/" if $outdir !~m!/$!;
my $tmppath = $outdir."tmp.json_$$";
open(TEMP,">$tmppath") or die "Couldn't write to a temporary file";
print TEMP "$json\n";
close TEMP;

my $outpath =$outdir.$outfile;
`mv $tmppath $outpath`;
