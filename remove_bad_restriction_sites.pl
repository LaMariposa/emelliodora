#!/usr/bin/perl

#remove_bad_restriction_sites.pl 
#script to remove read pairs if either one has bad restriction site, while allowing N reads for maintaining paired interleaved file
#M. Supple
#last modified 8 August 2016

#usage
#remove_bad_restriction_sites.pl <interleaved_fastq.gz> <proper_read_start>
	#<interleaved_fastq.gz> is a gzipped interleaved fastq
	#<proper_read_start> is expected start of each read (PstI=TGCAG)

#output
#gzipped interleaved fastq where read pairs with improper restriction site have been removed
#file is output in the working directory



use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my ($head1, $seq1, $qhead1, $qual1); 
my ($head2, $seq2, $qhead2, $qual2);

#die "$usage" unless (@ARGV == 2)

#read in command line arguments
my ($infile, $reseq)=@ARGV;
my $relen=length $reseq;

my @temp=split /\//, $infile;
my $outfile="clean_" . $temp[-1];
open(OUTFQ, "| gzip -c >$outfile");

#open input file
#open(INFQ, $infile) || die "can't open input file. $!\n";
open INFQ, "gunzip -c $infile |";

while ($head1=<INFQ>)
	{
		#read in information for both pairs
		$seq1=<INFQ>;
		$qhead1=<INFQ>;
		$qual1=<INFQ>;
		$head2=<INFQ>;
		$seq2=<INFQ>;
		$qhead2=<INFQ>;
		$qual2=<INFQ>;

		#check first sequence (keep if just N or has good restriction site
			#check if sequence is more than just single N
			if(length $seq1==2 || substr($seq1, 0, $relen) eq $reseq)
				{
					if(length $seq2==2 || substr($seq2, 0, $relen) eq $reseq)
						{
							print OUTFQ $head1;
                                                        print OUTFQ $seq1;
                                                        print OUTFQ $qhead1;
                                                        print OUTFQ $qual1;
                                                        print OUTFQ $head2;
                                                        print OUTFQ $seq2;
                                                        print OUTFQ $qhead2;
                                                        print OUTFQ $qual2;
						}
				}

	}

close INFQ;
close OUTFQ;
