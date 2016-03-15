#!/usr/bin/env perl 
#
# Manuel Rueda, PhD 09/2006
#
# Script for building the Hessian according to an elastic network model
# Physical Review Letters, 77, 9, 1905 (1996) Tirion, M.
#
# last 10/12/09

use strict;

if ( $#ARGV != 1 ) { print "Usage: $0 pdb output\n"; exit }

my $debug   = 0;
my $pdb     = $ARGV[0];      # input pdb file
my $out     = $ARGV[1];      # output hessian file
my $cutoff  = 8;             # cutoff for bonded atoms
my $cutoff2 = $cutoff**2;    # Computed once (avoids unncessary calculations)
my $kforce  = 0.25;          # In Kcal/mol.A2 units
my $kovacs  = 1;             # If Kovacs exponential function C=($dkvc/Rab)^6
                             # Proteins (2004) Volume 56, Issue 4, 661-668
my $cabond  = 1;             # J. Chem. Inf. Mod. (2009) 49 (3), pp 716–725
                             # In case we want to bond consecutive atoms (see below). Works with lineal/Kovacs
my $kbond = 350;             # Bond Constant for consecutive atoms
my $dkvc  = 3.8;             # Constant ~ distance between C-alphas
$cutoff = $cutoff2 = 999999
  if $kovacs;                #Just in case we forget manually to modificate the variable above

# Parse the input pdb
my ( @X, @Y, @Z, $acount, @chains );  # Note that arrays start at 1, not 0
open( PDB, "$pdb" ) or die "failed to open $pdb\n";
while (<PDB>) {
    next unless /^ATOM/ || /^HETAT/;
    $acount++;
    my ( $atom, $chain ) = readPDBLine($_);
    $X[$acount]      = $atom->{x};
    $Y[$acount]      = $atom->{y};
    $Z[$acount]      = $atom->{z};
    $chains[$acount] = $chain;
}
close PDB;

# Define N,N3
my $N  = $acount;    # Number of atoms
my $N3 = 3 * $N;     # Number of elements
print "CUT:$cutoff\n" if $debug;
print "N:$N\n"        if $debug;

# Hessian computation. Perl code adapted from Miyashita, O & Tama, F.
my ( $n, $m, $dx, $dy, $dz, $drdx, $drdy, $drdz, $r, $rr, @H, $k, $trace );

# Upper right matrix only (Sparse)
for $n ( 1 .. $N ) {
    for $m ( $n + 1 .. $N ) {

        $dx = $X[$n] - $X[$m];
        next if ( abs($dx) > $cutoff );
        $dy = $Y[$n] - $Y[$m];
        next if ( abs($dy) > $cutoff );
        $dz = $Z[$n] - $Z[$m];
        next if ( abs($dz) > $cutoff );
        $rr = $dx * $dx + $dy * $dy + $dz * $dz;
        next if ( $rr > $cutoff2 );

        $r = sqrt($rr);

        # The order of the below assignments DOES matter
        $k = $kforce if ( !$kovacs );
        $k = ( ( $dkvc / $r )**6 ) * $kforce if $kovacs;

       # From Wiki: "Tipically, two heavy atoms within 1.9 A are deemed to be covalently bonded"
       # We also want to keep connected atoms from gapped residues to avoid extreme fluctuations.
       # The assignement of K works as follows:
       # 'bonded atoms' and 'gapped i+1 consecutive (in PDB-file) atoms'
       #  (...and some 'collateral' -R end-atoms at d > 3.8 A from N

       # ** The method was initially 'optimized' for 'slightly' gapped proteins.
       # Later the web server was adapted to allow complete protein.
       # More recently, multi-chain objects (and heteroatoms) were also accepted.
       # We added the chain info (see above) to avoid that 2 consecutive atoms from different chains are Kbonded
       # they will still have the assigned elastic network k (kovacs/cutoff)

        $k = $kbond
          if ( $cabond
            && $chains[$n] eq $chains[$m]
            && ( $r <= 1.9 || ( $r > 3.8 && $m == $n + 1 ) ) );

        print "$n ($chains[$n]) $m ($chains[$m]) $k\n" if $debug;

        $drdx = $dx / $r;
        $drdy = $dy / $r;
        $drdz = $dz / $r;

        # Hii
        $H[ 3 * $n - 2 ][ 3 * $n - 2 ] += $k * $drdx * $drdx;
        $H[ 3 * $n - 2 ][ 3 * $n - 1 ] += $k * $drdx * $drdy;
        $H[ 3 * $n - 2 ][ 3 * $n ]     += $k * $drdx * $drdz;
        $H[ 3 * $n - 1 ][ 3 * $n - 1 ] += $k * $drdy * $drdy;
        $H[ 3 * $n - 1 ][ 3 * $n ]     += $k * $drdy * $drdz;
        $H[ 3 * $n ][ 3 * $n ]         += $k * $drdz * $drdz;

        # Hjj
        $H[ 3 * $m - 2 ][ 3 * $m - 2 ] += $k * $drdx * $drdx;
        $H[ 3 * $m - 2 ][ 3 * $m - 1 ] += $k * $drdx * $drdy;
        $H[ 3 * $m - 2 ][ 3 * $m ]     += $k * $drdx * $drdz;
        $H[ 3 * $m - 1 ][ 3 * $m - 1 ] += $k * $drdy * $drdy;
        $H[ 3 * $m - 1 ][ 3 * $m ]     += $k * $drdy * $drdz;
        $H[ 3 * $m ][ 3 * $m ]         += $k * $drdz * $drdz;

        # Hij
        $H[ 3 * $n - 2 ][ 3 * $m - 2 ] = -$k * $drdx * $drdx;
        $H[ 3 * $n - 2 ][ 3 * $m - 1 ] = -$k * $drdx * $drdy;
        $H[ 3 * $n - 2 ][ 3 * $m ]     = -$k * $drdx * $drdz;
        $H[ 3 * $n - 1 ][ 3 * $m - 2 ] = -$k * $drdy * $drdx;
        $H[ 3 * $n - 1 ][ 3 * $m - 1 ] = -$k * $drdy * $drdy;
        $H[ 3 * $n - 1 ][ 3 * $m ]     = -$k * $drdy * $drdz;
        $H[ 3 * $n ][ 3 * $m - 2 ]     = -$k * $drdz * $drdx;
        $H[ 3 * $n ][ 3 * $m - 1 ]     = -$k * $drdz * $drdy;
        $H[ 3 * $n ][ 3 * $m ]         = -$k * $drdz * $drdz;
    }
}

# Print the Hessian as "i j non-zero-ij-matrix-element" format
# Requires an external diagonalization program
# left-bottom corner matrix will be filled by external code
open( OUT, ">$out" ) or die "Cannot open $out\n";
for $n ( 1 .. $N3 ) {
    for $m ( $n .. $N3 ) {
        if ( defined $H[$n][$m] && $H[$n][$m] != 0 ) {
            printf( OUT "%6d %5d %15.12f\n", $n, $m, $H[$n][$m] );
        }
    }
}
close OUT;

sub readPDBLine {
    my ($line) = shift;
    my $newAt = {};
    $newAt->{name}      = substr( $line,          12, 4 );
    $newAt->{altloc}    = substr( $line,          16, 1 );
    $newAt->{residueId} = substr( $line,          17, 9 );
    $newAt->{x}         = substr( $line,          30, 8 );
    $newAt->{y}         = substr( $line,          38, 8 );
    $newAt->{z}         = substr( $line,          46, 8 );
    $newAt->{occ}       = substr( $line,          54, 6 );
    $newAt->{Bfact}     = substr( $line,          60, 6 );
    $newAt->{charge}    = $newAt->{occ};
    $newAt->{type}      = substr( $newAt->{name}, 1,  1 );

    # Extracting the chain info
    # Heteroatoms will belong to the same original chain
    my $chain = substr( $newAt->{residueId}, 4, 1 );
    return $newAt, $chain;
}
