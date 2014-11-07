#!/usr/bin/perl -w

package jsonReporter;

our $VERSION="1.0a";

use strict;
use warnings;
use DBI;
use LWP;
use JSON;
use XML::Simple;
use Data::Dumper;
use Time::Local;
use Time::Piece;
use POSIX qw(strftime);
use constant DEBUG=>0;

$| = 1;

my @Status = ("FAILED","KILLED","PREP","RUNNING","SUSPENDED","SUCCEEDED","OK");
my %topstats = (pending=>"pending",completed=>"completed",failed=>"failed"); # Binding internal definitions to external ones (which may change in a future)

use constant FAILED=>0;
use constant KILLED=>1;
use constant PREP=>2;
use constant RUNNING=>3;
use constant SUSPENDED=>4;
use constant SUCCEDED=>5;
use constant OK=>6;

use constant JSN=>"text/json";
use constant XML=>"text/xml";
use constant AGENT=>"SeqprodApp/0.1beta";

my $recentCutoff = 7 * 24 * 60 * 60; # week
my $recentCutoffTime = "week";
# query the seqware metadb
my $outfile;
my $outdir;

# OOZIE WEBSERVICE:
my $web_basedir = "11000/oozie/v1/job/";
my $web_bulkdir = "11000/oozie/v1/jobs?len=";
my $web_info    = "show=info";
my $web_defi    = "show=definition";

=head2 SUMMARY
 
 jsonReporter module is a placeholder for parameters used by seqprod_reporter web scripts which produce workflow run data in json format

=head2 USAGE
 
 To be updates soon
 
=cut


sub new {
    my $class = shift;
    my @args = @_;
    my $self = bless { dbhost    => $args[0],
                       dbname    => $args[1],
                       username  => $args[2],
                       password  => $args[3],
                       webservice=> $args[4],
                       devmode   => $args[5],
                       timemode  => $args[6]
                       }, ref $class || $class;    # will set finfiles flag to 1 if we have it (and if we do, genotype file and fingerprint popup will be created)
}

sub getSWData {
 my $self = shift;
 my $dsn = "DBI:Pg:dbname=$self->{dbname};host=$self->{dbhost}";
 my $dbh=DBI->connect($dsn, $self->{user}, $self->{password}, {RaiseError => 1});
 warn "connected to $self->{dbname}\n";

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
               wr.status_cmd
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

 my $records = $self->getWRData();
 print STDERR scalar(keys %{$records})." records received from webservice\n" if DEBUG;

 warn "Query complete, processing\n";
 my %results;
 my %statushash = ();
 while(my @row = $sth->fetchrow_array) {
        my ($sampleName, $workflowName, $workflowVersion, $statusCmd) = @row; 
        next if (!$sampleName || !$workflowName || !$workflowVersion || !$statusCmd);
        next if (!$records->{$statusCmd});

        my $crTime = parse_timestamp_oozie($records->{$statusCmd}->{ctime});
        my $lmTime = parse_timestamp_oozie($records->{$statusCmd}->{lmtime}) ;

        my $timeDiff = $currentTime - $crTime->[0];
	next if ($timeDiff > $recentCutoff);

        # Need to convert record's status in one of the three supported statuses here:
        my $currentStatus = $topstats{pending};
        if ($records->{$statusCmd}->{status} eq $Status[FAILED] || $records->{$statusCmd}->{status} eq $Status[KILLED]) {
          $currentStatus = $topstats{failed};
        } elsif ($records->{$statusCmd}->{status} eq $Status[OK] || $records->{$statusCmd}->{status} eq $Status[SUCCEDED]) {
          $currentStatus = $topstats{completed};
        }        


        $statushash{$records->{$statusCmd}->{status}}++;
        if (!$results{$currentStatus}){$results{$currentStatus}=[];}
        $statusCmd = (defined $statusCmd && $statusCmd=~/^\W*/) ? $' : ""; 
        push(@{$results{$currentStatus}},{sample    => $self->{devmode} ? $sampleName."[$records->{$statusCmd}->{status}]" : $sampleName,
                                          workflow  => $workflowName,
                                          version   => $workflowVersion,
                                          status    => $records->{$statusCmd}->{status}, #$status,
                                          #status_cmd=> $statusCmd, # will completely retire in a future
                	                  wrun_id   => $statusCmd, 
					  crtime    => $self->{timemode} ? $crTime->[0]*1000 : $crTime->[1],
					  lmtime    => $self->{timemode} ? $lmTime->[0]*1000 : $lmTime->[1]
                                          });
        
 }
 print STDERR "Seen these statuses:\n" if DEBUG;
 map{print STDERR "$_:\t$statushash{$_}\n"} (keys %statushash) if DEBUG;

 $sth->finish;
 $dbh->disconnect;

 return \%results;
}

#================================================================================================
# Subroutine for converting oozie time (i.e. Thu, 01 Jan 2009 02:00:00 GMT) to localtime (in sec)
#================================================================================================
sub parse_timestamp_oozie {
 my $tmpstmp = shift @_;
 if ($tmpstmp =~ /^(...),\s+(..)\s+(...)\s+(\d+)\s+(\d\d):(\d\d):(\d\d)\s+(\S+)/) {
                my $time = Time::Piece->strptime($tmpstmp, "%a, %d %b %Y %T %Z");
                # time format: 2014-09-23 14:44:10.032 (if timemode not set) 
                my $offset = $time->localtime->tzoffset;
                $time += $offset;
                print "Localized: ".$time."\n" if DEBUG;
                my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
                my $secs = timelocal($sec, $min, $hour, $mday, $mon, $year);

                # String should be formatted as: '2014-09-23 14:44:10.02'
                if ($year >= 1900){$year-=1900};
                my $now_string = strftime "%Y-%m-%e %T.00", $sec,$min,$hour,$mday,$mon,$year;
                $now_string =~s/\- /\-/;
                return [$secs,$now_string];
 }
 warn "Couldn't parse date [$tmpstmp]\n";
 return [0,"NA"];
}


#================================================
# Function that gets either total number of steps 
#    or steps done so far from a webservice
#================================================
sub countActions {
 
 my($self,$url,$mime) = @_;
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

#============================================================================
# Function for getting data from webservice, we get lmtime, crtime and status
#============================================================================
sub getWRData {
 my $self = shift;
 my $url = join(":",($self->{webservice},$web_bulkdir))."5"; # Initial request is to get the total number of entries
 my $req;

 print STDERR "Requesting $url\n" if DEBUG;
 $req = HTTP::Request->new(GET=>$url);
 $req->content_type(JSN);

 my $ua = LWP::UserAgent->new;
 $ua->agent(AGENT);
 my $res = $ua->request($req);

 # First, get the total number of records to bypass pagination
 my $total_records;
 if ($res->is_success) {
    print STDERR "JSON data received\n" if DEBUG;
    my $json_string = $res->content;
    my $j_container = decode_json $json_string;
    $total_records = $j_container->{total};
    $total_records ||= 0;
    print STDERR "Got ".$total_records." as total records from JSON\n" if DEBUG;
 } else {
   print STDERR $res->status_line, "\n"
 }

 if ($total_records == 0) {
    return undef;
 }

 # If we have more than 0 records, get the data
 $url = join(":",($self->{webservice},$web_bulkdir)).$total_records;
 print STDERR "Requesting $url\n" if DEBUG;
 $req = HTTP::Request->new(GET=>$url);
 $req->content_type(JSN);

 my $res_all = $ua->request($req);
 my $datachunks = {};

 if ($res_all->is_success) {
    print STDERR "JSON data received\n" if DEBUG;
    my $json_string_all = $res_all->content;
    my $j_container_all = decode_json $json_string_all;
    print STDERR "Got ".scalar(@{$j_container_all->{workflows}})." workflows from JSON\n" if DEBUG;
    foreach my $blah(@{$j_container_all->{workflows}}) {
        # time format: 2014-09-23 14:44:10.032 (if timemode not set)       
        $datachunks->{$blah->{id}} = {ctime => $blah->{createdTime},
                                      lmtime=> $blah->{lastModTime},
                                      status=> $blah->{status}};
    }
    print STDERR "Got ".scalar(keys %{$datachunks})." records from JSON\n" if DEBUG;    
 } else {
   print STDERR $res->status_line, "\n"
 }

 return $datachunks;
}

1;
