=b
This script helps you to convert a lmp data file to maxent fortran code input file (00input.dat).
use maxent_lurm.sh to submit your job. partition and other setting should be modifed for your job.
=cut

use strict;
use warnings;
use Data::Dumper;
use POSIX;
use Cwd;

open my $database ,"< HEC.data";#data file you want to apply maxent
my @data =<$database>;
close $database;
map { s/^\s+|\s+$//g; } @data;
my %para = (
    natom => "",
    ntype => "",
	xlo => "",
    xhi => "",
    ylo => "",
    yhi => "",
    zlo => "",
    zhi => "",
    coords => []
);
my $counter = 0;
for (@data){
    chomp;
    ####atoms###
    if(/(\d+)\s+atoms/){ 
        $para{natom} = $1;
    }
	####atom types###
    elsif(/(\d+)\s+atom types/){ 
        $para{ntype} = $1;
    }
    ####CELL_PARAMETERS###
    elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+xlo\s+xhi/){
        $para{xlo} = $1;
        $para{xhi} = $2;
    }
    elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+ylo\s+yhi/){
        $para{ylo} = $1;
        $para{yhi} = $2;
    }
    elsif(/([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+zlo\s+zhi/){
        $para{zlo} = $1;
        $para{zhi} = $2;
    }
    #1 1 4.458517505863 1.201338326940 0.873835074284
    elsif(/(\d+)\s+(\d+)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)\s+([+-]?\d*\.*\d*)$/){
        $counter++;
		chomp ($1,$2,$3,$4,$5);
		my $id = $1;
        my $type = $2;
		my $x = $3;
        my $y = $4;
        my $z = $5;
        my $temp = join(" ",(int($counter),$type, $x, $y, $z));
        #print "$temp\n";
        push @{$para{coords}},$temp;
    }
}#one data file    
my $coords = join("\n",@{$para{coords}});
chomp $coords;

my @box = (
"$para{xlo} $para{xhi}",
"$para{ylo} $para{yhi}",
"$para{zlo} $para{zhi}"
);
my $box = join("\n",@box);
chomp $box;

my $here_doc =<<"END_MESSAGE";
# The input file for MaxEnt developed by Prof. Shin-Pon Ju
####### begin below
#!the total MC iteration times,te total initial random search times
5000 10
#! PBC in x, y, and z
True True True 
#! rlist for cell list for the third nearest neighbour atoms
7.5
#!rdf peak(5) found by Ovito
4.5
6.5
7.5
8.6
9.5
#!total atom number and element type number for HEA
$para{natom} $para{ntype}
#!box information (xlo xhi...), only rectangular one is supported.
$box
#!coordinates
$coords
END_MESSAGE

open(FH, "> 00input.dat") or die $!;
print FH $here_doc;
close(FH);
