#!/usr/bin/env perl
#
# Manuel Rueda 09/2006 
#
# Script for Monte Carlo simulation in eigenvectors space
# Structure, Volume 15, Issue 5, 565-575, 16 May 2007
#
# Last modified 05/10/10

use strict;

if ( $#ARGV != 4 ) { print "$0 eval-force #frames #modes Temp seed\n"; exit }

my $fforce  = $ARGV[0];    # Eigenvector file
my $nframes = $ARGV[1];    # ~ Frames
my $itot    = $ARGV[2];    # Number of modes used
my $TEMP    = $ARGV[3];    # Temperature
my $seed    = $ARGV[4];    # Seed for the random

# $nmc is the number of MC movements
my $nmc = 100000;

# Initial seed
srand($seed);              #if we need to reproduce numbers

# Buffer of lines
my $buffer = int( 7900 * ( 5 / $nframes ) );

# fundamental constants
use constant BOLTZ => 0.00198717;

#
my $KbT  = BOLTZ * $TEMP;
my $BETA = 1 / $KbT;

# read force constants (Kcal/mol.A^2) ->@STIFF
my @STIFF = read_force($fforce);    #start at 0

# XCONF0 has the equilibrium parameters for each mode, set to 0.0
my @XCONF0 = ();
for my $i ( 0 .. ( $itot - 1 ) ) { $XCONF0[$i] = 0; }    # (although not necessary since empty=0, forced to be 0)

my @XCONFT = @XCONF0;                                    # Initial configuration of the try (at equilibrium)

# agr stands for the agressivity, the value
# of the maximum allowed displacement along the first mode. Adjust it to have
# around 40% acceptance. The other displacements are adjusted automatically
my $agr = 2.75 * sqrt( $KbT / $STIFF[0] );               #empirically calibrated

# scale the movements
my @SCAL = scaling();

# Energy for all modes at equilibrium = 0
my @ENER = ();
for my $i ( 0 .. ( $itot - 1 ) ) { $ENER[$i] = 0; }      # although not necessary since empty=0 (forced to be 0)

# ENER0 starting energy
my $ener0 = 0;

# Total energy at the begining
my $enertot = $ener0;

###################################################
#       Now starts MC Metropolis run              #
#       Select 1 mode randomly                    #
###################################################
#
my ( $dx, $dd, $tdx, $itt, $ittp, $delta, $enerbase, $enex, $etry, $newline, $iflag, $counter, $accepted );    #Logical value=0
for my $iter ( 1 .. $nmc ) {                                                                                   #Start LOOP
    # $itt from 0 to $itot-1
    # Problems if $itt==$itot (almost imposible, no cases although $nmc=1.000.000.000)

    # itt is the randomly selected mode to be changed
    $itt = int( rand($itot) );

    # note we sample both + and - displacements
    # $dx number between -1 and 1
    $dx = rand(2);
    $dx = $dx - 1;

    # we multiply it by the scaling factor
    $dx = $dx * $SCAL[$itt];

    # Energy of all the modes except the changed
    $enerbase = $enertot - $ENER[$itt];    # 0 in the first iteration

    # Displacement of the try, taking into account the previous for the same mode
    $tdx = $XCONFT[$itt] + $dx;

    # enex has the new energy of the step after changing itt variable
    $dd   = $tdx - $XCONF0[$itt];          # x-x0 -> $XCONF0[$it] shoukd be 0 in this case.
    $enex = $STIFF[$itt] * $dd * $dd;

    # This is the total energy of the try
    $etry = $enerbase + $enex;

    # Do metropolis test
    # Stochastic processes use a random sampling procedure to search conformational space
    # The score (energy) is calculated in each step and compared to the previous. If the new energy is lower
    # the step is accepted, otherwise the result is treated probabilistically by a Boltzmann mechanism
    # The higher the temperature the higher the likehood the step is accepted
    # PS: Another way to generate the DX foreach mode could be from 'normal random values'
    #    (see Box-Muller transform and Ziggurat algorithm) without any further Metropolis test.
    # PS2: The more modes we have the more time the system takes for equilibration since Dx in high modes will generate bigs $delta
    # that won't be accepted soonish. First moves will be low freq and progressively the others will be incorporated.
    # In 'normal-use' $buffer>10 & itot<500 system equilibrated before printing)
    # PS3: Better if avoiding subroutines/sending variables-references(slower code)
    $delta = $etry - $enertot;
    if ( $delta < 0 || exp( -$BETA * $delta ) > rand() ) {

        # accepted configuration
        $iflag = 1;
        $accepted++;
        $counter++;

        # New total energy
        $enertot = $etry;

        # New energy for the mode
        $ENER[$itt] = $enex;

        # New DX for the mode
        $XCONFT[$itt] = $tdx;

        # Print
        if ( $counter == $buffer ) {
            $ittp = $itt + 1;    #since arrays start at 0
            print "\*\*\*\*\n";
            print "ITERATION:$iter,MODE:$ittp,ETOT:$enertot,DX:$dx,CUM-DX($ittp):$tdx\n";
            print "SEED:$seed,NMODES:$itot\n";
            print "-------------------------------------------------------------------\n";
            $newline = 0;

            # DX Foreach mode
            for (@XCONFT) {
                $newline++;
                if ( $newline != 10 ) {
                    printf "%10.4f", $_;
                } else {
                    printf "%10.4f\n", $_;
                    $newline = 0;
                }
            }
            $counter = 0;
            print "\n\n" if ( $newline != 0 );
            print "\n"   if ( $newline == 0 );
        }

        #
    } else {
        $iflag = 0;
    }
}    # END LOOP

#
my $aratio = ( $accepted / $nmc ) * 100;
print "-------------------\n";
print "RATIO(%):$aratio\n";
print "-------------------\n";

############################################
# END MC
###########################################
sub scaling {
    my @SCAL;
    $SCAL[0] = $agr;
    for my $i ( 1 .. ( $itot - 1 ) ) {
        my $rat = $STIFF[0] / $STIFF[$i];
        $SCAL[$i] = sqrt( $agr * $agr * $rat );
    }
    return @SCAL;
}

sub read_force {
    my @FORCE;
    my $file = shift;
    my @LINE;
    my $counter;
    open( EVAL, "$file" ) || die "Cannot open $file\n" if $file !~ /\.gz$/;
    open( EVAL, "zcat $file |" ) || die "Cannot open $file\n" if $file =~ /\.gz$/;
    while (<EVAL>) {
        chomp;
        $counter = 1 if (/\*\*\*\*/);
        if ( $counter == 1 && /[0-9]/ ) {

            # just in case, keep only numbers
            @LINE = grep { /[0-9]/ } ( split( / +/, $_ ) );

            # Force constants divided by 2
            # E=0.5*K*x^2
            # print "$LINE[-1]\n";
            push @FORCE, $LINE[-1] / 2;
            $counter = 0;
        }
    }
    close EVAL;
    return @FORCE;    #starts at 0
}
