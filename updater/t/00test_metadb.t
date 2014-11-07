#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use Test::More 'no_plan';
use DBI;
#
# A simple test script that checks if we can connect to SeqWare DB
#
chdir $Bin;
use lib "$Bin/../blib";
use jsonReporter::ConfigData;

my $username = jsonReporter::ConfigData->config('username');
my $password = jsonReporter::ConfigData->config('password');
my $dbhost   = jsonReporter::ConfigData->config('dbhost');
my $dbname   = jsonReporter::ConfigData->config('dbname');



# We just need to see if we could connect to DB

my $result = &testDB;

if (!$result || $result ne "OK") {
 ok( 0, "No connection could be established with valid DB" );
} else {
 ok( 1, "All is good, we have a valid config" );
}

sub testDB {
 my $dsn = "DBI:Pg:dbname=$dbname;host=$dbhost";
 my $dbh=DBI->connect($dsn, $username, $password, {RaiseError => 1});
 return $dbh ? "OK" : "ERROR";
}

