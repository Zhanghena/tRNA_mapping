#!/usr/bin/env perl
#########################################################################
# File Name: removePrecursor.pl
# Author: Zhenghena
# mail: zhanghena@gmail.com
# Created Time: Thu 12 Jul 2018 09:42:00 AM CST
#########################################################################

use strict;
use warnings;
use Getopt::Long;

sub usage {
	print <<"END_USAGE";
Usage: perl $0
	--sam
	--bed
	--fq1
	--fq2
	--out1
	--out2
	--pre
	--mature
Example: perl $0 --sam *.sam --bed *.bed --fq1 *_R1.fq --fq2 *_R2.fq --out1 *_R1.fq --out2 *_R2.fq --preSam *.sam --matureSam *.sam
Only remove Precursor reads, not remove intron reads.
END_USAGE
	exit;
}

my ($sam,$bed,$fq1,$fq2,$out1,$out2,$pre,$mature);
GetOptions (
	'sam=s'      =>\$sam,
	'bed=s'      =>\$bed,
	'fq1=s'      =>\$fq1,
	'fq2=s'      =>\$fq2,
	'out1=s'     =>\$out1,
	'out2=s'     =>\$out2,
	'pre=s'      =>\$pre,
	'mature=s'   =>\$mature,
) or usage();
usage() if (!$sam or !$bed or !$fq1 or !$fq2 or !$out1 or !$out2 or !$pre or !$mature);

open SAM,    "< $sam"   or die "cannot open $sam\n";
open BED,    "< $bed"   or die "cannot open $bed\n";
open FQ1,    "< $fq1"   or die "cannot open $fq1\n";
open FQ2,    "< $fq2"   or die "cannot open $fq2\n";
open OUT1,   "> $out1"  or die "cannot open $out1\n";
open OUT2,   "> $out2"  or die "cannot open $out2\n";
open PRE,    "> $pre"   or die "cannot open $pre\n";
open MATURE, "> $out2"  or die "cannot open $mature\n";

my %tRNA  = ();
my %read  = ();
my @entry = ();

while(<BED>){
	chomp $_;

	my @id = split/\t/,$_;
	$tRNA{$id[3]} = $id[2]-$id[1];  #length
}

while(<SAM>){
	chomp $_;

	my @id = split/\t/,$_;
	my $start  = $id[3];
	my $len    = length($id[9]);
	my @cigar = split(/(I|M|X|D)/, $id[5]);
		my $I = 0;
		my $D = 0;
		for(my $i=0; $i < scalar(@cigar); $i++){
			if($cigar[$i] eq "I"){ $I += $cigar[$i-1];}
			if($cigar[$i] eq "D"){ $D += $cigar[$i-1];}
		}
	my $end    = $start + $len -1 + $D - $I;

	if($start > 100 && $end <= $tRNA{$id[2]}-96){
		$read{$id[0]} = 0;
		print PRE $_."\n";
	}
	else{
		print MATURE $_."\n";
	}
}

while(<FQ1>){
	chomp;
	push @entry, $_;
	if (scalar(@entry) == 4) {
		my ($id, $seq, $plusLine, $qual) = @entry;
		@entry = ();
		my @ID = split/(\@| )/,$id;
		if(exists $read{$ID[2]}){
			print OUT1 "$id\n$seq\n$plusLine\n$qual\n";
		}
	}
}


while(<FQ2>){
	chomp;
	push @entry, $_;
	if (scalar(@entry) == 4) {
		my ($id, $seq, $plusLine, $qual) = @entry;
		@entry = ();
		my @ID = split/(\@| )/,$id;
		if(exists $read{$ID[2]}){
			print OUT2 "$id\n$seq\n$plusLine\n$qual\n";
		}
	}
}
