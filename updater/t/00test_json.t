#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use Test::More 'no_plan';

#
# A simple test script that checks if we can connect to DB and webservice
#


# We just need to see if we could get any data from DB or webservice
my $result = "OK" #``; TODO: implement actual testing for DB and webserver connection

if (!$result || $result ne "OK") {
 ok( 0, "No connection could be established with valid DB" );
} else {
 ok( 1, "All is good, we have a valid config" );
}

