package NMA;

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Cwd;

=head1 NAME

    ALiBERO::NMA - Package for NMA calculation

=head1 SYNOPSIS

  use ALiBERO::NMA

=head1 DESCRIPTION


=head1 AUTHOR

Written by Manuel Rueda, PhD

=cut

=head1 METHODS

=cut

sub en_nma_mrc {
    my ($argsSub) = @_;

    my $N   = $main::args{ncpu};
    my $cmd = '';
    my %local_best_children =
      %::best_children;    # To avoid changing the original hash

    my $exe_dir  = '/pro/alibero/ALiBERO/bin';
    my $diaghess = "$exe_dir/diaghess_mrc_100_server"
      ;    #2500 atoms Compiled (crow) static 11/4/08 printing 100
    my $mc        = "$exe_dir/m-mc-eigen.pl";
    my $anim      = "$exe_dir/m-pca_anim_mc.pl";
    my $nma       = "$exe_dir/nmanuaa.pl";
    my $input_pdb = "best_child.pdb";
    my $f_mov     = 'mov';
    my $seed =
      int( rand(12379987497) )
      ; # srand only uses integer part. We need different $seed when elitism = on

    # We keep memory of the previous dir
    my $prev_dir = cwd();
    chdir( $argsSub->{receptor_dir} );

    # Creating global paths
    if ( $argsSub->{gen} > 1 ) {
        for my $element ( 1 .. $main::best_children{N} ) {
            $local_best_children{$element} =
              "$argsSub->{receptor_dir}/v$main::best_children{$element}.ob";
        }
    }

    # First we are going to copy the parents into MRC0
    # If Gen > 1 then docking results are copied via IO.pm
    mkdir('MRC0');
    my $tmpCount = 0;
    foreach my $element ( 1 .. $main::best_children{N} ) {
        $tmpCount++;
        $cmd =
"$main::rsync $local_best_children{$element} MRC0/RECEPTOR_$tmpCount.ob";
        system("$cmd") == 0 or die("failed to execute: $!\n");
    }

    # Defining here the number of MRC per children
    # We can have things like 10/3 or 25/4 so we add + 1 just in case
    $N = $N - $tmpCount;
    my $N_div_n_children = int( $N / $main::best_children{N} ) + 1;
    my $N_cum            = $tmpCount + 1;
    my $N_end            = $tmpCount + $N_div_n_children;

    # We skip the nma only if nParents = nChildren
    if ( $main::best_children{N} != $main::args{ncpu} ) {

        # EN-NMA (non-ICM version)
        foreach my $MRC ( 1 .. $main::best_children{N} ) {

            #We change the location into the RECEPTORs dir
            chdir( $argsSub->{receptor_dir} );
            mkdir("MRC$MRC");
            chdir("MRC$MRC");

            # We need to go from icm-object to non-H pdb for the NMA
            open( ICM, "| $main::icm $main::dev_null" );
            print ICM "read object \"$local_best_children{$MRC}\"\n";
            print ICM "write pdb a_1.//!h* \"$input_pdb\"\n";    # HA
                #print ICM "write pdb a_1.//ca* \"$input_pdb\"\n"; # CA
            print ICM "q\n";
            close ICM;

            # First Step NMA, build the Hessian Matrix
            # but only do it if the object is different, if not reuse
            print
"$main::prompt Building the Hessian for the child $main::best_children{$MRC}\n"
              if $main::args{debug} > 3;
            $cmd = "perl $nma $input_pdb hessian.dat $main::dev_null";
            system("$cmd") == 0 or die("failed to execute: $!\n");

            # Second Step, diagonalize Hessian
            print
"$main::prompt Diagonalizing the Hessian for the child $main::best_children{$MRC}\n"
              if $main::args{debug} > 3;
            $cmd = "$diaghess $main::dev_null";
            system("$cmd") == 0 or die("failed to execute: $!\n");
            system("/bin/gzip -f eigenvec.dat") == 0
              or die("failed to execute: $!\n");
            unlink <hessian.dat>;

            #Now generate cartesian structures
            print
"$main::prompt Creating $N_div_n_children cartesian MRCs (1-2 extra) for the child $main::best_children{$MRC}\n"
              if $main::args{debug} > 3;
            my $traj_dir = "TRAJ_${N_div_n_children}";
            mkdir $traj_dir;

# First we generate the projections for each mode
# Imp the space between $seed > $traj_dir/$f_mov
# Unfrequently, N gets N - 1 mrc. We add an extra 'ad-hoc' conformation just to make sure we get N total
            my $N_extra = $N_div_n_children + 1;
            $cmd =
"perl $mc eigenvec.dat.gz $N_extra $main::params{nmodes} $main::params{temperature} $seed > $traj_dir/$f_mov";
            system("$cmd") == 0 or die("failed to execute: $!\n");

# Now we use the projections to generate cartesian structures. Because of historic reasons, -pout is in AMBER CRD format
            $cmd =
"perl $anim -i $traj_dir/$f_mov -n $main::params{nmodes} -pout $traj_dir/$input_pdb -evec eigenvec.dat.gz -pdb $input_pdb $main::dev_null";
            system("$cmd") == 0 or die("failed to execute: $!\n");

            # We go from AMBER CRD to pdb format
            from_crd_to_pdb( $input_pdb, $traj_dir, "$input_pdb.x.gz",
                $N_div_n_children );

            # Tethering step with respect previous
            my $n_tmp = 1;
            for my $mrc ( $N_cum .. $N_end ) {
                $cmd =
"$main::params{macrodir}/Tether.icm $local_best_children{$MRC} $traj_dir/$input_pdb.$n_tmp $main::dev_null";
                system("$cmd") == 0 or die("failed to execute: $!\n");

                rename( "Tethered.ob", "RECEPTOR_${mrc}.ob" );
                print "$main::prompt OK => conformer $mrc tethered\n"
                  if $main::args{debug} > 3;
                $n_tmp++;
            }
            $N_cum = $N_cum + $N_div_n_children;
            $N_end = $N_end + $N_div_n_children;
        }
    }
    else {
        print "$main::prompt Skipping NMA step\n"
          if ( $main::args{debug} || $main::args{verbose} );

    }

    # We go back to the previous folder
    chdir $prev_dir;
}

sub from_crd_to_pdb {
    my $pdbref   = shift;
    my $traj_dir = shift;
    my $crd      = shift;
    $crd = "$traj_dir/$crd";
    my $N = shift;
    my ( @crds, @seqs );

    # Read Reference sequence from PDB
    open( PDBREF, "$pdbref" ) || die "cannot open $pdbref\n";
    while ( defined( my $line = <PDBREF> ) ) {
        next unless $line =~ /^ATOM/ || $line =~ /^HETATM/;
        chomp $line;
        my ($atom) = readPDBLine($line);
        push @seqs, $atom->{residueId};
    }
    my $natom  = scalar @seqs;
    my $natom3 = $natom * 3;

    # Read AMBER CRD file
    open( CRD, "/bin/zcat $crd |" ) || die "cannot open $crd\n";
    while ( defined( my $line = <CRD> ) ) {
        next if $line =~ /MC generated trajectory/;
        chomp $line;
        push @crds, grep { /[0-9]/ } ( split /\s+/, $line );
    }

    # Now print a pdb file foreach conformation
    my $a = 0;
    my $c = 0;
    for ( my $conf = 1 ; $conf <= $N ; $conf++ ) {
        my $b = 0;
        open( CONF, ">$traj_dir/$pdbref.$conf" )
          || die "cannot open $pdbref.$conf\n";
        for ( $a = $c ; $a < ( $c + $natom3 ) ; $a = $a + 3 ) {

            # if any difference with ptraj is @ -0.000 (ptraj 0.000)
            print CONF "$seqs[$b]";
            printf CONF "%8.3f%8.3f%8.3f %5.2f %5.2f\n", $crds[$a],
              $crds[ $a + 1 ], $crds[ $a + 2 ], 0, 0;
            $b++;
        }
        close CONF;
        $c = $a;
    }
}

sub readPDBLine {
    my ($line) = shift;
    my $newAt = {};
    $newAt->{residueId} = substr( $line, 0, 30 );
    return $newAt;
}
1;
