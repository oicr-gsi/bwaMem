#!/usr/bin/perl -w

use strict;
use FindBin '$Bin';

use lib "$Bin/install_util";
use File::Copy;
use jsonReportInstall;

  my $build = jsonReportInstall->new
    (module_name => 'jsonReporter',
     license     => 'perl',
     dist_author        => 'Peter Ruzanov <pruzanov@oicr.on.ca>',
     dist_abstract      => 'jsonReporter webscript for getting data from Seqware DB and Oozie webserver at OICR',
     requires => {
                  'perl'          => '5.14.0',
                  'CGI'     => '2.5',
                  'Config'  => '',
                  'JSON'       => '2.4',
                  'DBI' => '1.6',
                  'LWP' => '4.8',
                  'XML::Simple' => '2.0',
                  'Time::Local' => '1.1',
                  'Time::Piece' => '1.1',
                  'Getopt::Std' => '1.0',
                  'POSIX'       => '1.1'
                 }
    );

  $build->create_build_script;


print STDERR <<END;

**** WARNING *****************************************************************************************
This is rather unflexible setup, it requires user to specify a top-level directory for installation,
usually \$HOME/public_html or /var/www. The module and all the scripts with .html documents and images
will be installed in directory specified at the installation time. Changing ownership for the reports dir
will be required (apache needs to write into this directory)
ALSO, R needs to be installed! R scripts won't run if there's no Rscript in the PATH
******************************************************************************************************
END
