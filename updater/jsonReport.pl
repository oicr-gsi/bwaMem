#!/usr/bin/perl
use warnings;
use strict;

use DBI;
use JSON;
use LWP;
use XML::Simple;
use Data::Dumper;
use Time::Local;
use Getopt::Std;
use constant DEBUG=>0;
$| = 1;

# OOZIE WEBSERVICE:
#my $webservice  = "http://hsqwstage-node2.hpc.oicr.on.ca";
my $webservice  = "http://hsqwprod-node1.hpc.oicr.on.ca"; 
my $web_basedir = "11000/oozie/v1/job/";
my $web_info    = "show=info";
my $web_defi    = "show=definition";

my @Status = ("FAILED","KILLED","PREP","RUNNING","SUSPENDED","SUCCEEDED","OK");

use constant FAILED=>0;
use constant KILLED=>1;
use constant PREP=>2;
use constant RUNNING=>3;
use constant SUSPENDED=>4;
use constant SUCCEDED=>5;
use constant OK=>6;

use constant JSN=>"text/json";
use constant XML=>"text/xml";
use constant AGENT=>"SeqprodApp/0.1alpha";

my $recentCutoff = 7 * 24 * 60 * 60; # week
my $recentCutoffTime = "week";
# query the seqware metadb
my $username = '*****';
my $password = '*****';
my $dbhost = 'hsqwprod-db1.hpc.oicr.on.ca';
my $dbname = 'hsqwprod_seqware_meta_db';
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
    die "Please run as seqwareJsonReport.pl [week|month|year|decade] [output_dir]\n";
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

my $dsn = "DBI:Pg:dbname=$dbname;host=$dbhost";
my $dbh=DBI->connect($dsn, $username, $password, {RaiseError => 1});
warn "connected to $dbname\n";

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $currentTime = timelocal($sec,$min,$hour,$mday,$mon,$year);
my $query = "WITH RECURSIVE workflow_run_processings (workflow_run_id, processing_id) AS (
             SELECT wr.workflow_run_id, p.processing_id from workflow_run wr JOIN processing p
              ON (p.workflow_run_id is not null and wr.workflow_run_id = p.workflow_run_id)
              or (p.workflow_run_id is null and wr.workflow_run_id = p.ancestor_workflow_run_id)
             UNION
	     SELECT p.workflow_run_id, pr.parent_id FROM workflow_run_processings p JOIN processing_relationship pr
	      ON p.processing_id = pr.child_id), 
	      total_workflow_run_ius AS (
	     SELECT wr.workflow_run_id, i.ius_id FROM workflow_run wr 
	      JOIN ius_workflow_runs iwr ON wr.workflow_run_id = iwr.workflow_run_id 
	      JOIN ius i ON iwr.ius_id = i.ius_id
	     UNION 
	     SELECT wrp.workflow_run_id, i.ius_id FROM workflow_run_processings wrp 
	      JOIN processing_ius pi ON wrp.processing_id = pi.processing_id 
	      JOIN ius i ON pi.ius_id = i.ius_id)
             SELECT 
               wr.workflow_run_id,
               s.name AS sample_name,
               w.name,
               w.version,
               wr.status,
               wr.status_cmd,
               wr.create_tstmp,
               p.last_modified
              FROM workflow_run AS wr 
              JOIN workflow AS w ON wr.workflow_id = w.workflow_id
              JOIN total_workflow_run_ius AS twri ON twri.workflow_run_id = wr.workflow_run_id
              JOIN ius AS i ON i.ius_id = twri.ius_id
              JOIN sample AS s ON i.sample_id = s.sample_id
              LEFT OUTER JOIN (
                SELECT workflow_run_id, MAX(update_tstmp) AS last_modified
                FROM processing
                GROUP BY workflow_run_id ) AS p ON p.workflow_run_id = wr.workflow_run_id
              WHERE wr.status!='null'";
	        
my $sth = $dbh->prepare($query);
$sth->execute();

warn "Query complete, processing\n";
my %results;
my %statushash = ();
while(my @row = $sth->fetchrow_array) {
        my ($workflowRunID, $sampleName, $workflowName, $workflowVersion, $status, $statusCmd, $createTime, $lastmodTime) = @row;
               
        next if (!$workflowRunID || !$sampleName || !$workflowName || !$workflowVersion || !$createTime);
        $lastmodTime ||= $createTime;

        #my $ctime_value = scalar($currentTime);
        #my $ltime_value = parse_timestamp($lastmodTime);

        my $timeDiff = $currentTime - parse_timestamp($lastmodTime);
	next if ($timeDiff > $recentCutoff);
        if ($timemode) {
          $createTime = parse_timestamp($createTime) * 1000;
          $lastmodTime= parse_timestamp($lastmodTime) * 1000;
        }
        my $currentStatus = $status eq 'failed' || $status eq 'completed' ? $status : "pending";
        $statushash{$status}++;
        if (!$results{$currentStatus}){$results{$currentStatus}=[];}
        $statusCmd = (defined $statusCmd && $statusCmd=~/^\W*/) ? $' : ""; 
        push(@{$results{$currentStatus}},{sample    => $devmode ? $sampleName."[$status]" : $sampleName,
                                          workflow  => $workflowName,
                                          version   => $workflowVersion,
                                          status    => $status,
                                          status_cmd=> $statusCmd,
                	                  wrun_id   => $workflowRunID,
					  #accession => $workflowAccession,
					  crtime    => $createTime,
					  lmtime    => $lastmodTime});
        
}
print STDERR "Seen these statuses:\n" if DEBUG;
map{print STDERR "$_:\t$statushash{$_}\n"} (keys %statushash) if DEBUG;

$sth->finish;
$dbh->disconnect;

# Need to know progress
if (defined $results{pending} && scalar(@{$results{pending}}) > 0) {
 warn "Have pending workflows, will parse ".scalar(@{$results{pending}})." runs";
 RUN:
 foreach my $run (@{$results{pending}}) {
    if (!defined $run->{status} || $run->{status} ne 'running') {
     next RUN;
    }
          # FAKE STATUS IF NO PROGRESS AVAILABLE
          print STDERR "Using Fake status Data\n" if DEBUG;
          my $progress = "".int rand(100); 
          if (defined $run->{status_cmd} && $run->{status_cmd}=~/oozie/) { 
             # Create a request
             my $url = join(":",($webservice,$web_basedir)).$run->{status_cmd}."?";
             print STDERR "Getting oozie workflow run progress from $url\n" if DEBUG;
             my $steps = countActions($url.$web_defi,XML);
             my $done  = countActions($url.$web_info,JSN);
             $progress = int($done/$steps*100);
             print STDERR "Got real Progress value $progress\n"; 
          }
          $run->{progress} = $progress;
  } 
}
print STDERR "We have following types present in the results:\n" if DEBUG;
map{print STDERR $_."\n"} (keys %results) if DEBUG;

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

# Returns time in seconds
sub parse_timestamp {
 my $tmpstmp = shift @_;
 if ($tmpstmp =~ /^(....)-(..)-(..) (..):(..):(..)/) {
 print STDERR $tmpstmp."\n" if DEBUG;
                my $year = $1;
                my $mon = $2 - 1;
                my $mday = $3;
                my $hour = $4;
                my $min = $5;
                my $sec = $6;
                return timelocal($sec,$min,$hour,$mday,$mon,$year);
 } 
 warn "Couldn't parse date [$tmpstmp]\n";
 return 0;
}


#================================================
# Function that gets either total number of steps 
#    or steps done so far from a webservice
#================================================
sub countActions {

 my($url,$mime) = @_;
 my $req;

 if ($url && $mime) {
  print STDERR "Requesting $url\n" if DEBUG;
  $req = HTTP::Request->new(GET=>$url);
  $req->content_type($mime);
 } else {
  print STDERR "Cannot get content from the webservice\n";
  return undef;
 }

 my $ua = LWP::UserAgent->new;
 $ua->agent(AGENT);
 my $res = $ua->request($req);

 if ($mime eq JSN) {
   my $success = 0;
   if ($res->is_success) {
    print STDERR "JSON data received\n" if DEBUG;
    my $json_string = $res->content;
    my $j_container = decode_json $json_string;
    print STDERR "Got ".scalar(@{$j_container->{actions}})." records from JSON\n" if DEBUG;
    foreach my $blah(@{$j_container->{actions}}) {
      if ($blah->{externalId} && $blah->{externalId}=~/\d+/) {
         if ($blah->{status} eq $Status[OK] || $blah->{status} eq $Status[SUCCEDED]) {
           $success++;
         }
      }
    }
    } else {
      print STDERR $res->status_line, "\n"
    }
    return $success;
 } elsif ($mime eq XML) {
   my $steps = 1;
   if ($res->is_success) {
    print STDERR "XML data received\n" if DEBUG;
    my $xml_string = $res->content;
    my $x_container = XMLin($xml_string);
    $steps = scalar(keys %{$x_container->{action}}) - 1;
    $steps = 1 unless $steps >= 1;
   } else {
     print STDERR $res->status_line, "\n";
   }
   return $steps;
 }
 print STDERR "Format Unsupported\n";
 return 0;
}

