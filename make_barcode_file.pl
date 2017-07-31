#!/usr/bin/perl

#make_barcode_file.pl 
#script to generate a barcode file for axe demultiplexing from sample info spread sheets
#M. Supple
#last modified 22 April 2015

#usage
#make_barcode_file.pl <PlateName> <plate_db> <sample_db> <barcode_dir>

#output
#tab delimited file with information on axe demultiplexing


#improve by having it search the header to determine the column numbers
#ugly reformating of cell names

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my @plate_info;

#die "$usage" unless (@ARGV == 4)

#read in command line arguments
my ($plateName, $plateCSV, $sampleCSV, $barcodePlateDir)=@ARGV;

#open output file
my $out_file=$plateName . "_barcodes.tab";
open(OUTBC, ">$out_file");

#get info on barcode plate from plates.csv
open(PLATE, $plateCSV) || die "can't open file with plate information. $!\n";

#print some header info
print OUTBC ("#plate=", $plateName, "\n#plateDB=", $plateCSV, "\n#sampleDB=", $sampleCSV, "\n#\n");


#read in header and print 
my $line = <PLATE>;
print OUTBC ("#", $line);
while ($line = <PLATE>)
	{
	  chomp $line;
	  my @fields = split ",", $line;
	  if ($fields[0] eq $plateName)
		{
		  @plate_info=@fields;
		  #print header info to file
		  print OUTBC ("#", $line, "\n#\n");
		}
	}


#read barcode plate layout into hash
#determine barcode plate file
my $bc_file = $barcodePlateDir;
my $bc_plate=$plate_info[1];
my @bcs = split "Adapter", $bc_plate;
if ($bcs[0] eq "PstI"){$bc_file .= "/pst1_plate_"}
	else {print "Uh oh, something is wrong--either the BARCODE in the Plate.csv file is not in the expected format, or a different enzyme was used."}
my @temp = split "Plate", $bcs[1];
my $platenum;
if ($temp[1]<10){$platenum="0"}
$platenum .= $temp[1];
$bc_file .= $platenum . ".tab";
print OUTBC ("#barcode file=",$bc_file, "\n#\n");
#read barcodes into hash
my %barcodes;
open(BARCODES, $bc_file) || die "can't find barcode file. $! \n";
while($line = <BARCODES>)
	{
	  chomp $line;
	  my ($bc1, $bc2, $well) = split /\t/, $line;
	  my $comb_bc=$bc1 . "\t" . $bc2;
	  $barcodes{$well} = $comb_bc;
	}

#extract sample info from the samples.csv
my @samples;
#open sample file
open(SAMPLE, $sampleCSV) || die "can't open file with sample information. $!\n";
#read in each line of the file, looking for samples on the plate
while($line = <SAMPLE>)
	{
	  chomp $line;
	  my @fields = split ",", $line;
	  if ($fields[12] eq $plateName)
		{
		  #add extra info to format coodinates like files
		  my $coords=$platenum . substr $fields[13],0,1;
		  my $num=substr $fields[13], 1;
		  if ($num < 10){$coords .= "0"}
		  $coords .= $num;
		  #look up barcodes in hash and print to file
		  print OUTBC ($barcodes{$coords},"\t", $fields[9], "\n");
		}
	}












































