package Help;

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

=head1 NAME

    ALiBERO::Help - Help file for the ALiBERO script

=head1 SYNOPSIS

  use ALiBERO::Help

=head1 DESCRIPTION


=head1 AUTHOR

Written by Manuel Rueda, PhD

=cut

=head1 METHODS

=cut

sub usage {

    # http://www.gsp.com/cgi-bin/man.cgi?section=3&topic=Getopt::Long
    my $version = shift;
    my %args    = ();
    GetOptions(
        'debug=s'   => \$args{debug},        #string
        'v'         => \$args{version},      #flag
        'verbose'   => \$args{verbose},      #flag
        'help|h'    => \$args{help},         #flag
        'man'       => \$args{man},          #flag
        'n=i'       => \$args{ncpu},         #numeric (integer)
        'i|input=s' => \$args{configfile}    #string (-i as in AMBER MD package)
               #'O'         => \$args{overwritte}    #flag

    ) or pod2usage( -exitval => 0, -verbose => 1 );

    # Control check
    if ( $args{version} ) { print "$version\n"; exit 0 }
    pod2usage( -exitval => 0, -verbose => 2 ) if $args{man};
    pod2usage( -exitval => 0, -verbose => 1 ) if $args{help};
    pod2usage( -exitval => 0, -verbose => 1 )
      if ( !$args{ncpu} || !$args{configfile} );
    pod2usage(
        -exitval => 0,
        -verbose => 1,
        -message => 'Option --i requires a config_file'
    ) if ( !-s $args{configfile} );
    pod2usage(
        -exitval => 0,
        -verbose => 1,
        -message => 'Option --n requires a positive integer'
    ) if ( $args{ncpu} <= 0 );    # Must be positive integer
    $args{debug} = 0 if !$args{debug};
    return %args;                 # NonRef
}
1;
