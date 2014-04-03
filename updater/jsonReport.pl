#!/usr/bin/perl
use warnings;
use strict;

use DBI;
use JSON;
use Data::Dumper;
use Time::Local;
use constant DEBUG=>1;
my $STATUSTAG = "pegasus";

my $recentCutoff = 7 * 24 * 60 * 60; # week
my $recentCutoffTime = "week";
# query the seqware metadb
#my $username = 'seqware';
my $username = 'pruzanov';
#my $password = 'seqware7420';
my $password = 'metaPHP';
#my $dbhost = 'sqwprod-db1.hpc.oicr.on.ca';
#my $dbhost = '10.2.0.56';
#my $dbhost = '10.2.0.44';
my $dbhost = 'localhost';
my $dbname = 'meta_db';
my $wrkey  = '~/cron/flyking-rsync-key';
my $wruser = 'seqware';
my $wrhost = 'pipedev.hpc.oicr.on.ca';

if (exists $ARGV[0])
{
  if ($ARGV[0] eq "month")
  {
    $recentCutoffTime = $ARGV[0];
    $recentCutoff = 31 * 24 * 60 * 60;
  }
  elsif ($ARGV[0] eq "year")
  {
    $recentCutoffTime = $ARGV[0];
    $recentCutoff = 365 * 24 * 60 * 60;
  }
  elsif ($ARGV[0] eq "decade")
  {
    $recentCutoffTime = $ARGV[0];
    $recentCutoff = 10 * 365 * 24 * 60 * 60;
  }
  elsif ($ARGV[0] ne "week")
  {
    die "Please run as seqwareJsonReport.pl [week|month|year|decade] > output.html\n";
  }
  # defaults are fine
}


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

while(my @row = $sth->fetchrow_array) {
        my ($sampleName, $workflowName, $workflowVersion, $status, $statusCmd, $createTime, $lastmodTime) = @row;
               
        next if (!$sampleName || !$workflowName || !$workflowVersion || !$createTime);
        $lastmodTime ||= $createTime;
        

        #my $ctime_value = scalar($currentTime);
        #my $ltime_value = parse_timestamp($lastmodTime);

        my $timeDiff = $currentTime - parse_timestamp($lastmodTime);
	warn "$createTime vs $lastmodTime\n" if DEBUG;
        #warn "$ctime_value vs $ltime_value\n\n" if DEBUG;
        #print Dumper($createTime)."\n".Dumper($lastmodTime)."\n";
        #exit;
	next if ($timeDiff > $recentCutoff);
        
        my $currentStatus = $status eq 'failed' || $status eq 'completed' ? $status : "pending";
        
        if (!$results{$currentStatus}){$results{$currentStatus}=[];}
        $statusCmd = $' if $statusCmd=~/^\W*/; #'
        warn "Status command:".$statusCmd."\n" if DEBUG;
        push(@{$results{$currentStatus}},{sample    => $sampleName,
                                          workflow  => $workflowName,
                                          version   => $workflowVersion,
                                          status    => $status,
                                          stcommand => $statusCmd,
                	                  #id => $workflowRunID,
					  #accession => $workflowAccession,
					  #host => $wrhost,
					  crtime    => $createTime,
					  lmtime    => $lastmodTime});
        
}

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

    warn "About to check for ".$STATUSTAG if DEBUG;
    if (defined $run->{stcommand} && $run->{stcommand}=~/$STATUSTAG/) {
	  warn "Checking ".$run->{workflow}." with ".$run->{stcommand}."\n" if DEBUG;
	  my @statusOutput = split(/\n/, `ssh -i $wrkey $wruser\@$wrhost $run->{stcommand}`);
	  if (defined $statusOutput[5] && $statusOutput[5] =~ /.*\( (.*%) \).*/) {
	    $run->{progress} = $1;
	  } else {
            $run->{progress} = "NA";
          }
	 }
   }
}

my $json = JSON->new();
$json = encode_json \%results;
print "$json\n";


sub parse_timestamp {
 my $tmpstmp = shift @_;
 if ($tmpstmp =~ /^(....)-(..)-(..) (..):(..):(..)/) {
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
