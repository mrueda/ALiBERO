package PBS;

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

=head1 NAME

    ALiBERO::PBS - Package for submitting VS jobs to PBS

=head1 SYNOPSIS

  use ALiBERO::PBS

=head1 DESCRIPTION


=head1 AUTHOR

Written by Manuel Rueda, PhD

=cut

=head1 METHODS

=cut

our $local_host_name = `/bin/hostname -f`;
chomp $local_host_name;
our $qsub_sdsc = '/opt/torque/bin/qsub';
our $qsub      = 'export SGE_ROOT=/gridware/sge; /usr/local/bin/qsub';
our $qstat     = 'export SGE_ROOT=/gridware/sge; /usr/local/bin/qstat';
our $qdel      = 'export SGE_ROOT=/gridware/sge; /usr/local/bin/qdel';
our $c_lab     = 'c-ablab.ucsd.edu';
our $ssh_c     = "ssh $c_lab";

sub submit_dockscan {
    my ($argsSub) = @_;
    my $cmd       = '';
    my $status    = 0;

# Desktop / Workstation Mode
# Ad hoc solution to submit as many jobs as params{ncpuhost} (we avoided using Perl modules to simplify the installation)
    if ( $main::params{pbs} eq 'off' ) {

# ICMHOME/rundock -l 5. for thoroughness = 5
#        $cmd =
#"echo '#!/bin/bash\nexport ICMHOME=$main::ICMHOME\n/bin/hostname > info.txt\n\$ICMHOME/rundock -o TEST -f 1 -t $main::params{nligands} -a -l $main::params{thoroughness} $main::dev_null &' > dockScan.bash\nchmod +x dockScan.bash";

        $cmd =
"echo '#!/bin/bash\nexport ICMHOME=$main::ICMHOME\n/bin/hostname > info.txt\nmv TEST.dtb TEST.tab\n\$ICMHOME/rundock -o TEST -f 1 -t $main::params{nligands} -a -l $main::params{thoroughness} $main::dev_null &' > dockScan.bash\nchmod +x dockScan.bash";

        system("$cmd") == 0 or die("failed to execute: $!\n");
        my $running_jobs =
          `/bin/ps -ef | /bin/grep rundock | /bin/grep -cv grep`;
        chomp $running_jobs;
        while ( $running_jobs >= $main::params{ncpuhost} ) {
            print
"$main::prompt Waiting: Jobs running ($running_jobs) Total CPUs ($main::params{ncpuhost})\n"
              if ( $main::args{verbose} || $main::args{debug} );
            sleep(60);    # Wait 60 seconds before next check
            $running_jobs =
              `/bin/ps -ef | /bin/grep rundock | /bin/grep -cv grep`;
            chomp $running_jobs;
        }
        system("./dockScan.bash") == 0 or die("failed to execute: $!\n");
        $running_jobs = `/bin/ps -ef | /bin/grep rundock | /bin/grep -cv grep`;
        chomp $running_jobs;
        print
"$main::prompt Submitted: Jobs running ($running_jobs) Total CPUs ($main::params{ncpuhost})\n"
          if ( $main::args{verbose} || $main::args{debug} );
    }

    # Triton SDSC (Last Tested 2011. Needs to be updated)
    elsif ( $main::params{pbs} eq 'triton' ) {

        my $machine    = "m1rueda\@triton-login.sdsc.edu";
        my $ssh        = "ssh $machine";
        my $remote_dir = '/home/m1rueda/TEST_ALiBERO/';

        # Using a temporary name for the remote dir
        my $id_job = time . substr( "00000$$", -5 );            # Id for the job
        my $remote_full_dir = $remote_dir . "tmp_" . $id_job;

        my $tmp_bash_name = $argsSub->{repeat};
        $tmp_bash_name =~ s#\S+/(\S+)/(\S+)/(\S+)/(\S+)#$1_$2_$3_$4#;
        $tmp_bash_name =~ s#[AEIOU]##ig;
        $tmp_bash_name =
          'v' . $tmp_bash_name;    # because PBS may not like numbers at pos 1

        make_triton_pbs_cmd(
            {
                cmd_name   => $tmp_bash_name,
                n_ligands  => $main::params{nligands},
                local_dir  => $argsSub->{repeat},
                remote_dir => $remote_full_dir
            }
        );

        # It is IMP the / on final folder for the scp
        # scp -q does not print STDOUT
        # -B to run in batch, trying to avoid tty supended
        $cmd =
"scp -r -B -q $argsSub->{repeat}/ $machine:$remote_full_dir/ $main::dev_null";
        system("$cmd") == 0 or die("failed to execute: $!\n");
        $cmd =
"$ssh 'cd $remote_full_dir/ ; $qsub_sdsc $tmp_bash_name $main::dev_null'";
        system("$cmd") == 0 or die("failed to execute: $!\n");
    }

    # ABAGYAN's LAB CLUSTER.
    elsif ( $main::params{pbs} eq 'bluefish' ) {

        # Name for remote folder
        my $id_job = time . substr( "00000$$", -5 );    # Id for the job
        my $remote_full_dir = $main::params{cfsdir} . "/tmp_" . $id_job;

        # Using a temporary name for bash script
        my $tmp_bash_name = $argsSub->{repeat};
        $tmp_bash_name =~ s#\S+/(\S+)/(\S+)/(\S+)/(\S+)#$1_$2_$3_$4#;
        $tmp_bash_name =~ s#[AEIOU]##ig;
        $tmp_bash_name = 'v'
          . $tmp_bash_name
          . ".bash";    # because PBS may not like numbers at pos 1

        make_bluefish_pbs_cmd(
            {
                cmd_name   => $tmp_bash_name,
                n_ligands  => $main::params{nligands},
                local_dir  => $argsSub->{repeat},
                remote_dir => $remote_full_dir
            }
        );

        # Copying the files
        $cmd =
"$main::rsync $argsSub->{repeat}/* $ENV{USER}\@$c_lab:$remote_full_dir/ $main::dev_null";
        print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
        system("$cmd") == 0 or die("failed to execute: $!\n");

        # Submitting the jobs to PBS remotely
        # We need to capture the Id
        $cmd = "$ssh_c '$qsub $remote_full_dir/$tmp_bash_name'";
        print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
        my $status = `$cmd`;
        my ($job_id) = $status =~
          /Your job (\d+)/;    # important to have parenthesis at ($job_id)
        print "$main::prompt PBS Job Id ($job_id)\n" if $main::args{debug} > 3;

# every now and then 5-40% of the jobs fail because of /cfs not-mounted or not writtable
# ad-hoc solution => checking and re-submitting
        _check_Eqw_jobs_bluefish(
            {
                job_id     => $job_id,
                cmd_name   => $tmp_bash_name,
                remote_dir => $remote_full_dir,
            }
        );
    }
}

sub make_triton_pbs_cmd {

    # Note: This was tested many months ago and may not work as of today
    my ($argsSub) = @_;
    my $ICMHOME   = '/home/cxedwards/icm/icmd';
    my $log       = 'vls.icm.log';

    # Print PBS file
    open( PBS, ">$argsSub->{cmd_name}" )
      || die("Cannot open $argsSub->{cmd_name} pbs command file\n");
    print PBS "#\!/bin/bash\n";
    print PBS "#PBS -q small\n";
    print PBS "#PBS -l walltime=03:00:00\n";
    print PBS "export ICMHOME=$ICMHOME\n";

    #print PBS "cp \044PBS_O_WORKDIR/* \044PBSTMPDIR\n"; #TSRI
    #print PBS "cd \044PBSTMPDIR\n";                     #TSRI
    print PBS "cd $argsSub->{remote_dir}\n";
    print PBS "./MakeDock.icm $main::dev_null\n";
    print PBS
"\$ICMHOME/rundock TEST -f 1 -t $argsSub->{n_ligands} -o TEST -a -l $main::params{thoroughness} >& $log\n";

    # Now copying -f what is needed and discarding the rest #TSRI
    #print PBS "cp -f *.ou         \044PBS_O_WORKDIR\n"; #TSRI
    #print PBS "cp -f *answer*.ob  \044PBS_O_WORKDIR\n"; #TSRI
    #print PBS "cp -f $log         \044PBS_O_WORKDIR\n"; #TSRI
    #print PBS "rm -f *\n";                              #TSRI

    # And now from the master node we scp to Desktop
    #print PBS "cd \044PBS_O_WORKDIR\n";                #TSRI

    # We are going two steps in secure copy to make sure we got *answers* first
    # We put the / at the end destiny scp just in case
    print PBS
"scp -B -q *answer*.ob TEST.tab TEST_rec.ob $main::params{remoteuser}\@$local_host_name:$argsSub->{local_dir}/\n";
    print PBS
"scp -B -q *.ou  $log   $main::params{remoteuser}\@$local_host_name:$argsSub->{local_dir}/\n";

    # We erase the maps
    #print PBS "rm *map\n";
    print PBS "exit\n";
    close PBS;
    my $cmd = "chmod +x $argsSub->{cmd_name}";
    system("$cmd") == 0 or die("failed to execute: $!\n");
}

sub make_bluefish_pbs_cmd {

    my ($argsSub) = @_;
    my $log = 'vls.icm.log';

    # Print PBS file
    open( PBS, ">$argsSub->{cmd_name}" )
      || die("Cannot open $argsSub->{cmd_name} pbs command file\n");
    print PBS "#\!/bin/bash\n";
    print PBS
      "#\$ -S /bin/bash\n";    #Random behaviour because of non-w /home/user/
                               #print PBS "#\$ -j y -o $argsSub->{remote_dir}\n"
    ;
    print PBS "#\$ -j y -o /dev/null\n"
      ;    #We need /dev/null otherwise we can not delete dir
    print PBS "export ICMHOME=$main::ICMHOME\n";
    print PBS "cd $argsSub->{remote_dir}\n";
    print PBS "/bin/hostname > info.txt\n";
    print PBS "./MakeDock.icm $main::dev_null\n";

    print PBS "mv TEST.dtb TEST.tab\n"
      ;    # 02/13/13 Ad hoc solution to make it work with new ICM
    print PBS
"\$ICMHOME/rundock -o TEST -f 1 -t $argsSub->{n_ligands} -a -l $main::params{thoroughness} >& $log\n";
    print PBS
"$main::rsync *ou *answer*.ob TEST.tab TEST_rec.ob info.txt $log   $ENV{USER}\@$local_host_name:$argsSub->{local_dir}/\n"
      ;    # we need TEST_rec.ob in case we want to proccess individual results
    print PBS "rm *\n";                            #Careful with that
    print PBS "rmdir $argsSub->{remote_dir}\n";    #Careful with that
    print PBS "exit\n";
    close PBS;
    my $cmd = "chmod +x $argsSub->{cmd_name}";
    system("$cmd") == 0 or die("failed to execute: $!\n");
}

sub _check_Eqw_jobs_bluefish {

    # Note: Assuming that the jobs will enter immediatly to the queue
    my ($argsSub) = @_;
    my $status    = 0;
    my $cmd       = '';

    #my $sleep     = 20;
    my $sleep = 0;    # 0 for debugging

    #  We need at least 8 seconds to capture failed job
    sleep($sleep);

    # Check if the job is Eqw
    $cmd = "$ssh_c '$qstat -j $argsSub->{job_id} | /bin/grep -c error'";
    print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
    $status = `$cmd`;
    chomp($status);

    # STATUS = 0 is OK
    print "$main::prompt STATUS => $status\n" if $main::args{debug} > 3;

    # if the job was Eqw then delete and resubmit
    my $job_id = $argsSub->{job_id};
    until ( !$status ) {

        print "$main::prompt ** $job_id state is Eqw **\n"
          if $main::args{debug} > 3;

        # Delete
        print "$main::prompt Deleting $job_id\n" if $main::args{debug} > 3;
        $cmd = "$ssh_c '$qdel $job_id $main::dev_null'";
        system("$cmd") == 0 or die("failed to execute: $!\n");

        # Resubmit
        $cmd = "$ssh_c '$qsub $argsSub->{remote_dir}/$argsSub->{cmd_name}'";
        print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
        $status = `$cmd`;
        ($job_id) = $status =~
          /Your job (\d+)/;    # important to have parenthesis at ($job_id)
        print "$main::prompt PBS Job New Id ($job_id)\n"
          if $main::args{debug} > 3;

        # Wait X secs
        sleep($sleep);

        # Check again status
        print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
        $cmd    = "$ssh_c '$qstat -j $job_id | /bin/grep -c error'";
        $status = `$cmd`;
        chomp($status);
        print "$main::prompt NEW STATUS => $status\n" if $main::args{debug} > 3;
    }
}

sub copy_ligands_remote {

    # Name for remote folder
    my $remote_full_dir =
      $main::params{cfsdir} . "/tmp_LIGANDS_" . $main::params{id};

    # Creating an ICM file for remote index file
    ICM::create_index_file($remote_full_dir);

    # Copying the ligand files
    my $cmd =
"$main::rsync $main::params{sdf} $main::params{projdir}/LIGANDS/RemoteIndex.icm $ENV{USER}\@$c_lab:$remote_full_dir/ $main::dev_null";
    print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
    system("$cmd") == 0 or die("failed to execute: $!\n");
    $cmd = "$ssh_c $remote_full_dir/RemoteIndex.icm $main::dev_null";
    print "$main::prompt Executing $cmd\n" if $main::args{debug} > 3;
    system("$cmd");  # Ad-hoc solution to avoid communication issues with CentOs
                     #system("$cmd") == 0 or die("failed to execute: $!\n");
}
1;
