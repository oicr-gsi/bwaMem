#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use Test::More 'no_plan';
use LWP;
#
# A simple test script that checks if we can connect to Oozie webservice
#
chdir $Bin;
use lib "$Bin/../blib";
use jsonReporter::ConfigData;
use constant JSN=>"text/json";
use constant AGENT=>"SeqprodApp/0.1beta";

my $webservice = jsonReporter::ConfigData->config('webservice');
my $web_bulkdir = "11000/oozie/v1/jobs?len=5";

# We just need to see if we could get any data from Oozie webservice
my $url = join(":",($webservice,$web_bulkdir));
my $result = &testWS($url, JSN);

if (!$result || $result ne "OK") {
 ok( 0, "No connection could be established with Oozie websevice" );
} else {
 ok( 1, "All is good, we have a valid config" );
}

sub testWS {
 my($url,$mime) = @_;
 my $req;

 if ($url && $mime) {
  $req = HTTP::Request->new(GET=>$url);
  $req->content_type($mime);
 } else {
  print STDERR "Cannot get content from the webservice\n";
  return 0;
 }

 my $ua = LWP::UserAgent->new;
 $ua->agent(AGENT);
 my $res = $ua->request($req);
 my $success = "ERROR";
   if ($res->is_success) {
    warn "JSON data received";
    my $json_string = $res->content;
    $success = length($json_string) > 0 ? "OK" : "ERROR";
   }
 return $success;
}

