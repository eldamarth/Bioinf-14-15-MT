#!/usr/local/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config auto_help) ;
use File::Copy;
use Pod::Usage;

=head1 NAME

SNP Clustering

=head1 SYNOPSIS

perl ./SNPC.pl [options]

	Options:
		--rootdir	=>	Highest dir that contains all files. [C:/]
		--PHfile	=>	Polyhigh calls input data rootdir/./file. Substitute number by #.
		--SNPfile	=>	Core SNPs input data rootdir/./file. Substitute number by #.
		--groups    	=>  	Number of linkage groups to analyze.
		--outputdir	=>	Output root directory.	[rootdir/./LG/perc/]
		--percent	=> 	Similarity percentage cut	[97]
		--nind		=>	Number of individuals	[41]

=head1 DESCRIPTION

Perl software to cluster SNPs by chromosome based on their genotype.
Files for different chromosomes should be numbered.
Genotypes should be written in  genotype calls as 0,1,2 and -1 for No Call.
Number of individuals is set on 41, for different N, change $nind value.

=head1 AUTHOR

Modification by Miriam Paya of Chema Hidalgo's script - josemanuel.hidalgo@irta.cat

=cut
my $nind = "41";
my $PHfile = "";
my $SNPfile = "";
my $groups = "1";
my $rootdir = "C:/";
my $outputdir = "";
my $percent = "97";
my $PHchr;
my $SNPchr;
#------------Argument passing------------#
GetOptions (	"PHfile=s"   => \$PHfile,      # string
		"SNPfile=s"   => \$SNPfile,      # string
		"groups=i" => \$groups,    # numeric
		"rootdir=s"   => \$rootdir,      # string
		"outputdir=s"   => \$outputdir,      # string
		"nind=i" => \$nind,    # numeric
		"percent=i" => \$percent,    # numeric
	)
		or die ("Error in command line arguments\n");
				
#------------Creating Output Folder------------#
mkdir $rootdir.$outputdir or print STDERR "Couldn't create output folder $! \n";

for (my $group = 1; $group <= $groups; $group++) {
($PHchr = $rootdir.$PHfile) =~ s/\#/$group/;
($SNPchr = $rootdir.$SNPfile) =~ s/\#/$group/;

open PHFILEIN, $PHchr or die "Input file not found\n";
open SNPFILEIN, $SNPchr or die "Input file not found\n";

chdir $rootdir.$outputdir or print STDERR "Couldn't set output folder $! \n";
print "Creating Output folders for LG$group: \n";	# Modify folder name if desired.
my $outdir = "LG".$group;
mkdir $outdir or print STDERR "Couldn't create output folder $! \n";
chdir $outdir or print STDERR "Couldn't set output folder $! \n";
$outdir = $percent;
mkdir $outdir or print STDERR "Couldn't create output folder $! \n";
chdir $outdir or print STDERR "Couldn't set output folder $! \n";
open FILEOUT, ">", "unclusterout.txt" or die "Could not create output file";
open CLUSTEREDOUT, ">", "clusterout.txt" or die "Could not create output file";

#------------Script------------#
my $i;
my $a;
my $b;
my $het;
my $hom;
my $gaps;
my $invers;
my $perc;
my $percinv;
my @cluster;

my @PHarray = <PHFILEIN>;
close (PHFILEIN);
my @snparray = <SNPFILEIN>;
close (PHFILEIN);

printf FILEOUT $snparray[0];
splice(@snparray, 0, 1);

for ($i = 0; $i < scalar @PHarray; $i++) {
chomp $PHarray[$i] ;
$PHarray[$i] =~ tr/\r//d;  #Remove carriage returns ?
}
for ($i = 0; $i < scalar @snparray; $i++) {
chomp $snparray[$i] ;
$snparray[$i] =~ tr/\r//d; 
}

for ($i = 0; $i < scalar @snparray; $i++) {
	my $line = $snparray[$i]; 
	my @frase1 = split('\t', $line);
	my $id1 = $frase1[0];
	my $snpid1 = $frase1[0]."/".$frase1[1];
	push(@cluster, $line ); 
	for ($b = 0; $b < scalar @PHarray; $b++) {
		my @frase2 = split('\t', $PHarray[$b]);
		my $id2 = $frase2[0];
		if ( $id2 eq $id1 ) {
			next;
		}
		$het = 0;
		$hom = 0;
		$gaps = 0;
		$perc = 0;
		$invers = 0;
		for ($a = 3; $a < scalar @frase1; $a++)  # Check $a values match call positions
		{
			if (( $frase1[$a] eq "1" ) && ( $frase2[$a-2] eq "1" ))
			{
				$het++;
			}
			elsif (( $frase1[$a] eq $frase2[$a-2] ))
			{
				$hom++;
			}
			elsif  (( $frase1[$a] eq "0" ) && ( $frase2[$a-2] eq "2" ) || ( $frase1[$a] eq "2" ) && ( $frase2[$a-2] eq "0" ))
				{
					$invers++;
				}
			elsif (( $frase1[$a] ne "-1" ) && ( $frase2[$a-2] eq "-1" ) || ( $frase2[$a-2] ne "-1" ) && ( $frase1[$a] eq "-1" ))
				{
				$gaps++;
				}
			else
			{
			}
		}
		$perc=(($hom+$het+$gaps)/$nind)*100;
		$percinv=(($invers+$het+$gaps)/$nind)*100;
		if (( $perc > "$percent" ) || ( $perc eq "100" ))  {
			my $call = join("\t",$snpid1,"NORMAL",$PHarray[$b]);
			push(@cluster, $call); # Save similar SNPs in cluster array.
			if ( $snparray[$i] =~ $id1 ) {
				splice(@snparray, $i, 1) ;	# Remove clustered SNP from array
				$i--;
			}
		}
		elsif (( $percinv > "$percent" ) || ( $percinv eq "100" ))  
		{
			my $call = join("\t",$snpid1,"INVERS",$PHarray[$b]);
			push(@cluster, $call);	# Save similar SNPs in cluster array.
			if ( $snparray[$i] =~ $id1 ) {
				splice(@snparray, $i, 1) ;	# Remove clustered SNP from array
				$i--;
			}
		}
	}
	if ($cluster[-1] eq $line) {
		pop(@cluster)  # Removing unclustered snps
	}
}

foreach (@snparray) {
	 print FILEOUT "$_\n";
	 };

foreach (@cluster) {
	 print CLUSTEREDOUT "$_\n";
	 };

close FILEOUT;
close CLUSTEREDOUT;
}
