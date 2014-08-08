#!/usr/bin/perl -w

use warnings;
use CGI qw(:standard -debug escape);
use JSON;
use constant DEBUG=>0;

my $fileDir = "/.mounts/labs/PDE/web/html/spbreporter/";
my %files   = (week   => "week.json",
               month  => "month.json",
               year   => "year.json",
               decade => "decade.json");


if (my $q = param('range')) {
    chomp($q);
    my $file = $q.".json";
    warn "Got file $file";
    if ($files{$q}) {
       warn "Will try to load the file" if DEBUG;
       my $path = $fileDir.$files{$q};
        if (!-e $path) {
          $path = $fileDir.$files{week}; # If there's no requested file, return the week.json
        } 
        open(JSON,"<$path") or die "Couldn't read file [$path]";
        my $jsonBuffer;
        while (<JSON>){
         chomp;
         $jsonBuffer.=$_;
        }
        close JSON;
        &printjson($jsonBuffer);
    }
} else {
 warn "Parameter not detected" if DEBUG;
}

sub printjson {
 my $jsonText = shift @_;
 print header('application/json');
 print $jsonText;
}

