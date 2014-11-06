#!/usr/bin/perl -w

use strict;
use Cwd;
use FindBin '$Bin';

my $ROOTDIR = "spbreporter";
my $origdir = cwd;
my $homedir = "$Bin/..";

chdir $homedir or die "couldn't cd to $homedir: $!\n";
my $basedir   = shift @ARGV;
my $scriptdir = shift @ARGV;

$scriptdir ||= $basedir; # If we don't have scriptdir defined, use basedir

print STDERR "Starting installing into $basedir\n";
chop($basedir) if $basedir =~m!/$!;

# Make directory ntree
`mkdir -p $basedir/html/$ROOTDIR`;
`mkdir -p $basedir/cgi-bin/$ROOTDIR`;
`mkdir -p $scriptdir/$ROOTDIR/jsonReporter`;

# Copy all files
my $wd = `pwd`;
chomp($wd);
print STDERR "Working dir is $wd\n";

&copy_dir("$wd/cgi-bin",$basedir.'/cgi-bin/'.$ROOTDIR.'/',".pl"); 
&copy_dir("$wd/scripts",$scriptdir.'/'.$ROOTDIR.'/',".pl");
&copy_dir("$wd/blib/lib",$scriptdir.'/'.$ROOTDIR.'/',".pm");
&copy_dir("$wd/blib/lib/jsonReporter",$scriptdir.'/'.$ROOTDIR.'/jsonReporter/',".pm");

print STDERR <<END

***********************************************************************
Don't forget to do set up cronjobs for jsonReport.pl script (crontab -e) :

*/5 * * * * /u/pruzanov/seqprod_reporter/jsonReport.pl week /.mounts/labs/PDE/web/html/spbreporter/ 2> /dev/null
0 1 * * 0 /u/pruzanov/seqprod_reporter/jsonReport.pl month /.mounts/labs/PDE/web/html/spbreporter/ 2> /dev/null

***********************************************************************
END
;

# ===========================
# copy a dir content
# ===========================
sub copy_dir {
 my($src,$dest,$wildcard) = @_;
 if (-d $src && -d $dest) {
   opendir(DIR,$src) or die "Couldn't read from directory [$src]";
   my @files = $wildcard ? grep{/$wildcard$/} readdir DIR : grep{!/\.$|\.\.$/} readdir DIR;
   close DIR;

   foreach my $file (@files) {
     `cp $src/$file $dest/`;
   }
 } else {
   warn "Something wrong with the directories, cannot copy from $src";
 }
}
