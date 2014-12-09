#!/usr/bin/perl -w
#This program is used to uncompress the compressed hetfa and mask files.  
#by Mengyao Zhao, 2014-12-05

use strict;
use warnings;

my $sample = "uncompress";

if (@ARGV != 2) { 
	print ("Usage: uncompress.pl <reference.fa> <sample.ccomp.fa.gz>\n");
	exit;
}

if ($ARGV[1] =~ /(\w+\/)?(\S+)\.ccomp\.fa\.gz/) {
	$sample = $2;
} else {
	warn "Your uncompressed files will be uncompress.fa, uncompress.fai, uncompress.filter.fa and uncompress.filter.fai\n";
}

my $file = &gwhich ("uncomp") || die;

if (defined $1) {
	system("$file $ARGV[0] $ARGV[1] > $1$sample.fa");
} else {
	system("$file $ARGV[0] $ARGV[1] > $sample.fa");
}
system("samtools faidx $sample.fa");
system("gunzip -f $sample.ccompmask.fa.gz");
system ("mv $sample.ssompmask.fa $sample.filter.fa");
system("samtools faidx $sample.filter.fa");

sub which
{
	my $file = shift;
	my $path = (@_)? shift : $ENV{PATH};
	return if (!defined($path));
	foreach my $x (split(":", $path)) {
		$x =~ s/\/$//;
		return "$x/$file" if (-x "$x/$file");
	}
	return;
}

sub gwhich {
	my $progname = shift;
	my $addtional_path = shift if (@_);
	my $dirname = &dirname($0);
	my $tmp;
	chomp($dirname);
	if ($progname =~ /^\// && (-x $progname)) {
		return $progname;
	} elsif (defined($addtional_path) && ($tmp = &which($progname, $addtional_path))) {
		return $tmp;
	} elsif (-x "./$progname") {
		return "./$progname";
	} elsif (defined($dirname) && (-x "$dirname/$progname")) {
		return "$dirname/$progname";
	} elsif (($tmp = &which($progname))) {
		return $tmp;
	} else {
		return;
	}
}

sub dirname {
	my $prog = shift;
	return '.' if (!($prog =~ /\//));
	$prog =~ s/\/[^\s\/]$//g;
	return $prog;
}


