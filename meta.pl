#!/usr/bin/perl
# Author: Jason Liu
# Version: 1.0

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;
use File::Basename;
use File::Spec;

&main;
exit;

sub main
{
	my $help = 0;
	my %opts = ('help' => \$help);

	GetOptions(\%opts, 'snptest', 'method=i', 'dir=s', 'threshold:f', 'lambda:s', 'size:i', 'output=s', 'help');  
	pod2usage(1) if $help;
	## If no arguments were given, then allow STDIN to be used only
	## if it's not connected to a terminal (otherwise print usage)
 	#pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

	my $message;
	unless(defined $opts{method})
	{
    		$message = "Error! Method not sepecified.\n";
    		pod2usage(-exitval => 2, -msg => $message);
	}

	unless(defined $opts{dir})
	{
    		$message = "Error! Directory of cohorts not sepecified.\n";
    		pod2usage(-exitval => 2, -msg => $message);
	}

	unless(defined $opts{output})
	{
    		$message = "Error! Output file  not sepecified.\n";
    		pod2usage(-exitval => 2, -msg => $message);
	}

	unless(defined $opts{threshold})
	{
		$opts{threshold} = 0.5;
	}

	opendir DIR, $opts{dir} || die "Can't open directory:$!";
	my @files = sort(readdir DIR);
	my @cohorts;
	foreach my $file(@files) 
	{
		unless($file =~ /^\./)
		{
			push(@cohorts, "$opts{dir}$file");
		}
	} 
	
	my $com = "--method $opts{method} --cohort @cohorts";
	if(defined $opts{snptest})
	{
		$com = " --snptest ".$com;
	}
	
	if(defined $opts{threshold})
	{
		$com .= " --threshold $opts{threshold}"
	}

	if(defined $opts{lambda} || defined $opts{size})
	{
		my $info;
		$info = $opts{lambda} if(defined $opts{lambda});
		$info = $opts{size} if(defined $opts{size});
		
		my @lambdas;
		my @sizes;
		open INFO, $info || die "Can't find the cohorts information: $!";
                my $header = <INFO>;
		while(my $line = <INFO>)
		{
			chomp($line);
			my @item = split /\s/, $line;
			push(@lambdas, $item[1]) if($opts{lambda});
			push(@sizes, $item[2]) if($opts{size});
		}
		close INFO;

		$com .= " --lambda @lambdas" if($opts{lambda});
		$com .= " --size @sizes" if($opts{size});
	}
	
	$com = "./meta-stat-1.3 $com --output $opts{output}";
	print "$com\n";
	system($com);

	close DIR;
}


__END__

=head1 NAME

meta.pl -- A perl script to wrap the META.

=head1 SYNOPSIS

./meta.pl [options] [file ...]
 
=head1 DESCRIPTION

The meta.pl is used to wrap the META, which will ease the use of the META.

=head1 OPTIONS
 
=over 8
 
=item B<--snptest>

Optional, used if the output of SNPTEST is used.

=item B<--method>

Compulsory, 

method = 1, inverse-variance method for fixed-effects model;
method = 2, inverse-variance method for random-effects model;
method = 3, z-score combination method for fixed effects model.

=item B<--dir>

Compulsory, the path where the cohorts are stored.

=item B<--lambda>

Optional, the genomic control lambda for each cohort.

=item B<--size>

Optional, sample size for each cohort.

=item B<--output>

Compulsory, specify the name of output file.

=item B<--help>

Print a brief help message and exits.
 
=back
 
=head1 AUTHOR

Jason Liu <lt>jsliu@stats.ox.ac.uk<gt>.

=head1 COPYRIGHT

Copyright 2009, by Jason Liu at University of Oxford

=cut
