package ICM;

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use lib '/pro/alibero/ALiBERO/PBS';

=head1 NAME

    ALiBERO::ICM - Package for ICM tasks

=head1 SYNOPSIS

  use ALiBERO::ICM

=head1 DESCRIPTION

NB: With the exception of Tether.icm, we are not passing arguments to icm scripts. We are replacing strings inside

=head1 AUTHOR

Written by Manuel Rueda, PhD

=cut

=head1 METHODS

=cut

sub read_icb {
    my %best_children  = ();
    my $work_dir       = $main::params{projdir} . "/GEN_1";
    my $rec_dir        = $work_dir . "/RECEPTORS";
    my $icm_file       = $rec_dir . "/readICB.icm";
    my $order_for_dock = 'order4Dock.txt';

    # Preparate all the initial I/O (only Gen 1)
    mkdir("$rec_dir");

    # Create cmd file for ICM
    open( ICM, ">$icm_file" )
      || die("cannot open $icm_file\n");
    print ICM "#!$main::icm\n";
    print ICM "call _startup\n";
    print ICM "openFile \"$main::params{inputicb}\"\n";
    print ICM "set directory \"$rec_dir\"\n";
    print ICM "write column Name(Obj(a_*.)) \"$order_for_dock\"\n";
    print ICM "for i=1,Nof(a_*.)\n";
    print ICM " write object a_\$i.\n";
    print ICM "endfor\n";
    print ICM "quit\n";
    close ICM;
    my $cmd = "chmod +x $icm_file; $icm_file $main::dev_null";
    system("$cmd");

    # Number of best_children (order irrelevant here)
    chomp( my $n_best_children = `/bin/ls -1 $rec_dir/*.ob | wc -l` );
    $best_children{N} = $n_best_children
      ;    # Remember that we will only print a maximum in the log file

    my $element = 0;

    # The objects will be loaded according to ICM file
    foreach my $object (
        `/bin/sed '1,2d' $rec_dir/$order_for_dock | awk '{print \$1}'`)
    {      # It does not matter that the ICM markers go beyond n_best_children
        $element++;
        chomp $object;
        $best_children{$element} =
          $rec_dir . '/' . $object . '.ob';    # We'll need the full path later
    }
    return %best_children;                     #NonRef
}

sub make_dock {

    my ($argsSub) = @_;
    my $cmd = '';

    # Only if #1
    if ( $argsSub->{n_repeat} == 1 ) {

        my $makeMaps_file = 'MakeDock.icm';
        my $f_makemap     = "$main::params{macrodir}/MakeDock.icm";
        my $tmp_icm       = `/bin/cat $f_makemap`;

        # Substitutions
        my $lig_path = $main::params{inx};
        $lig_path =~ s/.inx//;
        $tmp_icm  =~ s/LIGAND_FILE/$lig_path/
          if $main::params{pbs} eq 'off';    # only once in file
        my $tmp_sbs =
            $main::params{cfsdir}
          . "/tmp_LIGANDS_"
          . $main::params{id} . "/"
          . $main::params{lig};
        $tmp_icm =~ s/LIGAND_FILE/$tmp_sbs/
          if $main::params{pbs} ne 'off';    # only once in file
        $tmp_icm =~ s/XXXXX/$argsSub->{rec_name}/;    # only once in file

        open( DOCKSCAN, ">./$makeMaps_file" )
          || die("can not open $makeMaps_file\n");
        print DOCKSCAN $tmp_icm;
        close DOCKSCAN;

        $cmd = "chmod +x $makeMaps_file";
        system("$cmd");
        $cmd = "./$makeMaps_file $main::dev_null";
        system("$cmd") if $main::params{pbs} eq 'off';
    }
    else {
        system("$main::rsync ../1/* .") == 0
          or die("failed to execute: $!\n");

        # sometimes the calcs already started and we get .BAK
        my @goners = glob "v*.bash *answers* *.ou info.txt *.BAK";
        foreach my $trash (@goners) {
            unlink $trash;
        }

    }
}

sub check_end_vls {
    my ($argsSub) = @_;
    my $ok        = 0;
    my $output    = "TEST_1-$main::params{nligands}.ou";
    my $answers   = 'TEST_answers1.ob';

    # The file answers must be in the dir sice is copied before .ou
    if ( -e "$argsSub->{repeat}/$output" && -e "$argsSub->{repeat}/$answers" ) {
        chomp( my $n_lines =
              `/bin/grep -m 1 -F -c FINISHED $argsSub->{repeat}/$output` );
        $ok = 1 if $n_lines == 1;
        $ok = 0 if $n_lines != 1;
    }
    return $ok;    #NonREF
}

sub hitlist {
    my ($argsSub)   = @_;
    my $cmd         = '';
    my $function    = uc( $main::params{function} );
    my $hit_file    = "$argsSub->{gen_dir}/Hitlist.icm";
    my $f_hitlist   = "$main::params{macrodir}/Hitlist.icm";
    my $hitlist_log = 'Hitlist.log';
    my $tmp_icm     = `/bin/cat $f_hitlist`;
    if ( $argsSub->{refin} eq 'off' )
    {              # we can not use $params{refinement} here
        ;
        $tmp_icm =~ s/XXXX/$argsSub->{gen}/;    #only once
    }
    if ( $argsSub->{refin} eq 'on' ) { # we can not use $params{refinement} here

# Note that we are including all the DOCK_{?,??,???} dirs instead
# of using the best + refined.
# If laziness > 0 it may happen that the Gen_Ref_X.icb has a different winner (rare) than Gen_X.icb.
        $hit_file = "$argsSub->{gen_dir}/Hitlist_Ref.icm";
        $tmp_icm =~ s/XXXX/Ref_$argsSub->{gen}/;    #only once
        $hitlist_log = 'Hitlist_Ref.log';
    }

    my $tables_kept;
    my %best_children = ();

    # Keeping a maximum of 10 Tables to avoid ICM complaints about #arrays
    if ( $main::args{ncpu} > 5 ) {
        $tables_kept = 5;
    }
    else {
        $tables_kept = $main::args{ncpu};
    }

    # We check if the #ligands is complete
    $tmp_icm =~ s/NLIGANDS/$main::params{nligands}/g;    #appears twice

    # Change the name of the Ligand inx and sdf file at Hitlist.icm
    my $lig_path = $main::params{inx};
    $lig_path =~ s/.inx//;
    $tmp_icm  =~ s/LIGAND_FILE/$lig_path/;               # only once in file
    $tmp_icm  =~ s/INX_NAME/$main::params{lig}/;         # only once in file

    # Replace the fitness function
    $function = 'NSA' if $function eq 'NSAPLUS';    # Ad-hoc solution (Note uc)
    $tmp_icm =~ s/FUNCTION/$function/;              # Appears once

    # Replace the distance_switch
    $tmp_icm =~ s/DISTANCE_CHECK/$main::params{rdistance}/;    # Appears once

    # Replace the # MRC
    $tmp_icm =~ s/NMRC/$main::params{mrc}/;                    # Appears once

    # Print Hitlist file
    open( HITLIST, ">$hit_file" ) || die("can not open $hit_file\n");
    print HITLIST $tmp_icm;
    close HITLIST;

    $cmd = "chmod +x $hit_file; $hit_file >& $hitlist_log";
    print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
    system("$cmd");

    # Parsing the values of best src/mrc combination
    $best_children{auc} =
`/bin/grep -m 1 -F -A2 'T_Results_Auc_Marker'       $hitlist_log | tail -1`;
    $best_children{nsa} =
`/bin/grep -m 1 -F -A2 'T_Results_Nsa_Marker'       $hitlist_log | tail -1`;
    $best_children{nsaplus} = $best_children{nsa};
    $best_children{score} =
`/bin/grep -m 1 -F -A2 'T_Results_SCORE_AV_Marker'  $hitlist_log | tail -1`;
    $best_children{con} =
      `/bin/grep -m 1 -F -A2 'T_Results_Con_Marker'  $hitlist_log | tail -1`;

    # Best best_children for mrc (order does not matter)
    chomp( my $n_best_children = `/bin/ls -1 v*.ob |wc -l` );
    $best_children{N} = $n_best_children;

    foreach my $element ( 1 .. $n_best_children )
    {    # It does not matter that the ICM markers go beyond n_best_children
        $best_children{$element} =
`/bin/grep -m 1 -F -A2 'T_Results_Rec_Name_Marker_${element}_END'  $hitlist_log | tail -1`;
        $best_children{$element} = substr( $best_children{$element}, 5, 9999 );
    }

    # We must initialize the rest to avoid warnings when printing
    foreach my $element ( $n_best_children + 1 .. $::max_children ) {
        $best_children{$element} = 0;
    }

    # Remove \n in Hash
    chomp %best_children;

    # recent error found when using Rdistance (X-windows connection)
    chomp( my $hitlist_error =
          `/bin/grep -m 1 -c "Can\'t connect to xDisplay" $hitlist_log` );

    #foreach  my $key ( keys %best_children ) {
    #print $key, " => #", $best_children{$key}, "#\n";
    #}

    return %best_children;    #NonRef
}

sub refinement {
    my ($argsSub) = @_;
    my $ref_file  = "$argsSub->{gen_dir}/Refine_Hitlist.icm";
    my $f_ref     = "$main::params{macrodir}/Refine_Hitlist.icm";
    my $ref_log   = 'Refine_Hitlist.log';
    my $tmp_icm   = `/bin/cat $f_ref`;

    #
    $tmp_icm =~ s/LIGAND_FILE/$main::params{inx}/;    # only once in file
    $tmp_icm =~ s/INX_NAME/$main::params{lig}/;       # only once in file
    $tmp_icm =~ s/NTOP/$main::params{ntop}/g;

    #
    open( REFN, ">$ref_file" ) || die("can not open $ref_file\n");
    print REFN $tmp_icm;
    close REFN;

    #
    my $cmd = "chmod +x $ref_file; $ref_file >& $ref_log";
    system("$cmd");
    chomp( my $n_ref_OK = `/bin/ls -1 rComplex*.ob | wc -l` );
    return $n_ref_OK;                                 #NonRef
}

sub create_index_file {

    # Name for remote folder
    my $remote_full_dir = shift;
    my $icm_file        = '';
    my $cmd             = '';
    $icm_file = "$main::params{projdir}/LIGANDS/LocalIndex.icm";
    $icm_file = "$main::params{projdir}/LIGANDS/RemoteIndex.icm"
      if $remote_full_dir;

    # Create cmd file for ICM
    open( ICM, ">$icm_file" )
      || die("cannot open $icm_file\n");
    print ICM "#!$main::icm\n";
    print ICM "call _startup\n";
    print ICM
"makeIndexChemDb \"$main::params{sdf}\" \"$main::params{inx}\" \"mol\" { \"ID\" } \n"
      if !$remote_full_dir;
    print ICM
"makeIndexChemDb \"$remote_full_dir/$main::params{lig}.sdf\" \"$remote_full_dir/$main::params{lig}.inx\" \"mol\" { \"ID\" } \n"
      if $remote_full_dir;
    print ICM "quit\n";
    close ICM;
    $cmd = "chmod +x $icm_file; $icm_file $main::dev_null" if !$remote_full_dir;
    system("$cmd") if !$remote_full_dir;
    $cmd = "chmod +x $icm_file" if $remote_full_dir;
    system("$cmd") if $remote_full_dir;
}
1;
