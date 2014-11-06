package jsonReportInstall;

use base 'Module::Build';
use strict;
use warnings;
use ExtUtils::MakeMaker 'prompt';

my %MY_PROPS = (basedir       => 'Directory for cgi-bin and html files for spbreporter ?',
                user          => 'user name for connecting to SeqWare DB',
                password      => 'password for SeqWare DB',
                dbhost        => 'Host for SeqWare DB',
                dbname        => 'SeqWare DB name',
                webservice    => 'URL of SeqWare webservice');

sub ACTION_build {
    my $self = shift;
    $self->depends_on('config');
    $self->SUPER::ACTION_build;
}

sub ACTION_reconfig {
    my $self = shift;
    $self->config_done(0);
    unless (Module::Build->y_n("Reuse previous configuration as defaults?",'y')) {
        for (keys %{$self->my_props}) {
            $self->config_data($_=>undef);
        }
    }
    $self->depends_on('config_data');
    warn "\n**Path reconfigured. Running \"Build clean\".\n";
    $self->ACTION_clean;
}

sub ACTION_config {
    my $self  = shift;
    local $^W = 0;

    my $prefix = $self->install_base || $self->prefix || '';

    return if $self->config_done;

    print STDERR "\n**Beginning interactive configuration**\n";

    my $props = $self->my_props;
    my %opts  = map {
        $_=>$self->config_data($_)
    } keys %$props;

    my $dire_warning = 0;
    my @keys = (keys %MY_PROPS);
    while (@keys) {
        my $key = shift @keys;
        my $conf_dir = $props->{$key} =~ /directory/i;

        $opts{$key} = prompt($props->{$key},$opts{$key});
        if ($conf_dir) {
            my ($volume,$dir) = File::Spec->splitdir($opts{$key});
            my $top_level     = File::Spec->catfile($volume,$dir);
        
            unless (-d $top_level) {
                next if Module::Build->y_n("The directory $top_level does not exist. Use anyway?",'n');
                redo;
            }

        }
        $opts{basedir} ||= "/var/www/";
     }
    

    for my $key (keys %opts) {
        $self->config_data($key=>$opts{$key});
    }

    $self->config_done(1);

    print STDERR "\n**Interactive configuration done. Run './Build reconfig' to reconfigure**\n";
        
}

sub ACTION_install {
    my $self = shift;
    my $prefix = $self->install_base || $self->prefix || '';

    $self->depends_on('config_data');
    $self->install_path->{basedir}    ||= $self->config_data('basedir');
    system("perl","install_util/install_files.pl",$self->install_path->{basedir});

    print STDERR "\n***INSTALLATION COMPLETE***\n";
}


sub config_done {
    my $self = shift;
    my $done = $self->config_data('config_done');
    $self->config_data(config_done=>shift) if @_;
    warn "NOTE: Run ./Build reconfig to change existing configuration.\n" if $done;
    return $done;
}

sub ACTION_config_data {
    my $self = shift;
    $self->depends_on('config');
    $self->SUPER::ACTION_config_data;
}


sub my_props {
    return \%MY_PROPS;
}

sub asString { return 'SwapReporter installer' }

1;
