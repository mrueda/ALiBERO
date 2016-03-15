package IO;

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Cwd qw(abs_path);
use Sys::Hostname;

=head1 NAME

    ALiBERO::IO - Package for IO subroutines

=head1 SYNOPSIS

  use ALiBERO::IO

=head1 DESCRIPTION


=head1 AUTHOR

Written by Manuel Rueda, PhD

=cut

=head1 METHODS

=cut

sub io_preparation {
    my ($argsSub) = @_;
    my $work_dir  = $main::params{projdir} . "/GEN_" . $argsSub->{gen};
    my $prev_dir  = $main::params{projdir} . "/GEN_" . $argsSub->{bgsf}
      if $argsSub->{gen} > 1;    #prev dir only needed if Gen > 1

    if ( $argsSub->{gen} > 1 ) {
        mkdir("$work_dir/RECEPTORS");
        for my $element ( 1 .. $main::best_children{'N'} ) {
            print
"$main::prompt Executing $main::rsync $prev_dir/v$main::best_children{$element}.ob $work_dir/RECEPTORS\n"
              if $main::args{debug} > 3;
            system(
"$main::rsync $prev_dir/v$main::best_children{$element}.ob $work_dir/RECEPTORS"
            );
        }
    }
}

sub io_copy_parents {
    my ($argsSub) = @_;
    my $status    = 0;
    my $work_dir  = $main::params{projdir} . "/GEN_" . $argsSub->{gen};
    my $prev_dir  = $main::params{projdir} . "/GEN_" . $argsSub->{bgsf}
      if $argsSub->{gen} > 1;

    # Copying the VS folder and renaming
    for my $element ( 1 .. $main::best_children{'N'} ) {

        print
"$main::prompt Executing cp -r $prev_dir/DOCK_$main::best_children{$element} $work_dir/DOCK_$element\n"
          if $main::args{debug} > 3; # cp instead of rsycn to avoid overwritting
        $status = system(
"cp -r $prev_dir/DOCK_$main::best_children{$element} $work_dir/DOCK_$element"
        );
        if ( ( $status >>= 8 ) != 0 ) {
            die
"$main::prompt We could not cp -r $prev_dir/DOCK_$main::best_children{$element} $work_dir/DOCK_$element\n";
        }
        for my $repeat ( 1 .. $main::params{repeat} ) {
            print
"$main::prompt Renaming $work_dir/DOCK_$element/${repeat}/RECEPTOR_$main::best_children{$element}.ob to $work_dir/DOCK_$element/${repeat}/RECEPTOR_$element.ob\n"
              if $main::args{debug} > 3;
            rename(
"$work_dir/DOCK_$element/${repeat}/RECEPTOR_$main::best_children{$element}.ob",
                "$work_dir/DOCK_$element/${repeat}/RECEPTOR_$element.ob"
            );
        }
    }
}

sub read_config_file {
    my $config_file = shift;
    my $ncpuhost    = `/bin/grep -c '^processor' /proc/cpuinfo`;
    chomp($ncpuhost);

    # These are the default values
    my %params = (
        pbs          => 'off',
        inputicb     => undef,
        projdir      => 'ALiBERO_TEST',
        cfsdir       => undef,
        temperature  => 300,
        function     => 'nsa',
        iauc         => 0,
        insa         => 0,
        iscore       => 0,
        icon         => 0,
        maxauc       => 100,
        maxnsa       => 100,
        maxscore     => -99,
        maxcon       => 1,
        maxgen       => 10,
        mrc          => 3,
        rdistance    => 'off',
        nligands     => 30,
        sdf          => undef,
        macrodir     => undef,
        refinement   => 'off',
        ntop         => 3,
        nmodes       => 100,
        repeat       => 1,
        laziness     => 0,                # 5 for cluster mode
        elitism      => 'on',
        thoroughness => '1.0',
        ncpuhost     => $ncpuhost,
        remoteuser   => $ENV{USER}
    );
    open( CONFIG, $config_file )
      || die "Cannot open $config_file config file\n";
    while ( defined( my $line = <CONFIG> ) ) {

        next if $line =~ /^\s*#/;
        chomp $line;
        $line =~ s/#.*//;                 # no comments
        $line =~ s/^\s+//;                # No leading white
        $line =~ s/\s+$//;                # No trailing white

        my @l_params = split /\s+/, $line;    # Overwritten with $line
                                              # Now we simplify naming
        my ( $var, $value ) = ( lc( $l_params[0] ), $l_params[1] );

        # Check user typos in parameters name
        my $param_syntak_ok = grep { $_ eq $var } keys %params;
        die "$main::error_prompt Parameter '$var' does not exist (typo?)\n"
          if !$param_syntak_ok;

        $value = 'off' if lc($value) eq 'no';     # For consistency
        $value = 'on'  if lc($value) eq 'yes';    # For consistency

        # Value assignation
        if (   $var eq 'inputicb'
            || $var eq 'projdir'
            || $var eq 'sdf'
            || $var eq 'macrodir' )
        {
            $params{$var} = abs_path($value);
        }
        elsif ( $var eq 'ncpuhost' ) {

            # forcing integer
            $params{$var} = abs($value);

        }
        elsif (
               $var eq 'pbs'
            || $var eq 'refinement'
            || $var eq 'function'
            || $var eq 'elitism'

          )
        {
            $params{$var} = lc($value);
        }
        elsif ( $var eq 'thoroughness' ) {
            $params{$var} = $value * 1.0;
        }

        else { $params{$var} = $value; }

    }
    close(CONFIG);

    # Below are a few internal paramaters that do not have default values
    my @i_names = split /\//, $params{sdf};
    $i_names[-1] =~ s/.sdf//;
    $params{lig}      = $i_names[-1];
    $params{inx}      = $params{projdir} . "/LIGANDS/" . $params{lig} . ".inx";
    $params{id}       = time . substr( "00000$$", -5 );
    $params{hostname} = hostname;
    $params{user}     = $ENV{USER};

    # Check if the important files exist and user typos in file location
    die "$main::error_prompt Check that the input files exist (or typos in $config_file)\n"
      if ( !-e $params{inputicb}
        || !-d $params{macrodir}
        || !-s $params{sdf} );

    die
"$main::error_prompt ICM does not allow Index file names containing - (minus) character\n"
      if ( $params{inx} =~ m/-/ );    # _Hitlist will fail if not
                                      # SDF
    die "$main::error_prompt sdf file does not have extension .sdf\n"
      if ( $params{sdf} !~ m/.sdf/ )
      ;    # To avoid the typical error on MakeDock.icm s/LIGAND_FILE
    die
"$main::error_prompt ICM does not allow SDF file names containing - (minus) character\n"
      if ( $params{sdf} =~ m/-/ );    # Hitlist will fail if not
    die
"$main::error_prompt Are you sure that you want to create a folder named $params{projdir} ?\n"
      if ( $params{projdir} =~ m/\.in$/ );

    # ProjDir
    die "$main::error_prompt Project dir $params{projdir} exists\n"
      if ( -d $params{projdir} );

    # Check that rundock exists
    die "$main::error_prompt Can't find ICM's $::ICMHOME/rundock\n"
      if ( !-f "$::ICMHOME/rundock" );

    # Miscellanea
    die
"$main::error_prompt Not allowed value for parameter pbs: $params{pbs} [off|triton|bluefish]\n"
      unless ( $params{pbs} eq 'off'
        || $params{pbs} eq 'triton'
        || $params{pbs} eq 'bluefish' );
    die
"$main::error_prompt You have not defined cfsdir parameter. Try something like /cfs/$ENV{USER}/\n"
      if ( $params{pbs} eq 'bluefish' && !$params{cfsdir} );
    die "$main::error_prompt rdistance only can have 'on' or 'off' values\n"
      unless ( $params{rdistance} eq 'on' || $params{rdistance} eq 'off' );
    die "$main::error_prompt ncpuhost must be lower or equal to $ncpuhost\n"
      if ( $params{ncpuhost} > $ncpuhost );

    # Initialize values for the fitness function
    if ( $params{function} eq 'auc' ) {
        $params{maxthreshold}  = $params{maxauc};
        $params{initthreshold} = $params{iauc};
    }
    elsif ( $params{function} eq 'nsa' || $params{function} eq 'nsaplus' ) {
        $params{maxthreshold}  = $params{maxnsa};
        $params{initthreshold} = $params{insa};
    }
    elsif ( $params{function} eq 'con' ) {
        $params{maxthreshold}  = $params{maxcon};
        $params{initthreshold} = $params{icon};
    }
    elsif ( $params{function} eq 'score' ) {

        # Using 'positive' values for the energy scores to simplify algorithm
        $params{maxscore}      = $params{maxscore} * -1;
        $params{iscore}        = $params{iscore} * -1;
        $params{maxthreshold}  = $params{maxscore};
        $params{initthreshold} = $params{iscore};
    }
    else {
        die
"$main::error_prompt Not allowed value for parameter function: $params{function} [auc|nsa|nsaplus|score|con]\n"

    }
    return %params;    #NonRef
}
1;
