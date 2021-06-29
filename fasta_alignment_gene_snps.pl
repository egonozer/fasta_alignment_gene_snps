#!/usr/bin/perl

use strict;
use warnings;

my $usage = "

fasta_alignment_gene_snps.pl [options] <alignment.fasta>

Identifies SNP positions in a mutiple alignment and notes whether they are
synonymous or non-synonymous. Also identifies frameshifts relative to the
reference (FS).

Only works in alignments against a single complete gene. Align using MUSCLE
or MAFFT or similar multiple aligner.

Options:
  -r    ID of the reference sequence
        (default: first sequence is assumed to be reference)

";

use Getopt::Std;
our ($opt_r);
getopts('r:');

die $usage unless @ARGV;

my $refid   = $opt_r if $opt_r;

my @seqs;
my $refseq;
open (my $in, "<$ARGV[0]") or die "ERROR: Can't open $ARGV[0]: $!\n";
my ($id, $seq);
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            $refid = $id if !$refid;
            if ($id eq $refid){
                $refseq = $seq;
            } else {
                push @seqs, ([$id, $seq]);
            }
        }
        $seq = "";
        $id = substr($line, 1);
        next;
    }
    $line =~ s/\s//g;
    $seq .= uc($line);
}
close ($in);
if ($id){
    if ($id eq $refid){
        $refseq = $seq;
    } else {
        push @seqs, ([$id, $seq]);
    }
    $seq = "";
} else {
    die "ERROR: No sequences in $ARGV[0]\n";
}
die "ERROR: reference ID $refid not found in the file\n" unless $refseq;

my %trans = translate();


foreach my $slice (@seqs){
    my ($qid, $qseq) = @{$slice};
    print "$qid";
    my ($rcod, $qcod) = ("") x 2;
    my ($rpos, $qpos) = (0) x 2;
    my $is_fs;
    for my $i (1 .. length($refseq)){
        my ($rbase, $qbase) = (substr($refseq, $i-1, 1), substr($qseq, $i-1, 1));
        if ($rbase ne "-"){
            $rpos++;
            $rcod .= $rbase;        
        }
        if ($qbase ne "-"){
            $qpos++;
            $qcod .= $qbase;
        }
        if (length($rcod) == 3 and length($qcod) == 3){
            my $end_fs;
            if ($is_fs){
                
                my $trcod = substr($refseq, $i-3, 3);
                my $tqcod = substr($qseq, $i-3, 3);      
                my $trpos = $rpos + 1;
                if ($tqcod =~ m/-/){
                    for my $j (reverse 0 .. 2){
                        if (substr($tqcod, $j, 1) ne "-"){    
                            $trpos -- ;
                        } else {
                            last;
                        }
                    }
                } else {
                    for my $j (reverse 0 .. 2){
                        if (substr($trcod, $j, 1) ne "-"){
                            $trpos --;
                        } else {
                            last;
                        }
                    }
                }
                print "\tENDFS$trpos";
                $is_fs =  0;
                $end_fs = 1;
            }
            my @rbases = split(//, $rcod);
            my @qbases = split(//, $qcod);
            my @out;
            my $basecount = $rpos - 3;
            for my $j (0 .. $#rbases){
                $basecount++;
                if ($rbases[$j] ne $qbases[$j]){
                    push @out, "$rbases[$j]$basecount$qbases[$j]";
                }
            }
            if (@out){
                print "\t", join(",", @out) unless $end_fs;
                my ($raa, $qaa) = ($trans{$rcod}, $trans{$qcod});
                if ($raa and $qaa){
                    if ($raa ne $qaa){
                        my $apos = $rpos/3;
                        print ",[$raa$apos$qaa]";
                    }
                }
            }
            ($rcod, $qcod) = ("") x 2;
            next;
        } elsif (length($rcod) == 3 or length($qcod) == 3) {
            $rcod = "" if length($rcod) == 3;
            $qcod = "" if length($qcod) == 3;
            next if ($is_fs);
            $is_fs = 1;
            
            my $trcod = substr($refseq, $i-3, 3);
            my $tqcod = substr($qseq, $i-3, 3);
            my $trpos = $rpos - 3;
            if ($tqcod =~ m/-/){
                for my $j (0 .. 2){
                    if (substr($tqcod, $j, 1) eq "-"){
                        $trpos += $j;
                        last;
                    }
                }
            } else {
                for my $j (0 .. 2){
                    if (substr($trcod, $j, 1) eq "-"){
                        $trpos += $j;
                        last;
                    }
                }
            }
            print "\tFS$trpos";
        }
    }
    if (length($rcod) == 3 and length($qcod) == 3){
        my $end_fs;
        if ($is_fs){  
            my $trcod = substr($refseq, length($refseq)-3, 3);
            my $tqcod = substr($qseq, length($refseq)-3, 3);      
            my $trpos = $rpos + 1;
            if ($tqcod =~ m/-/){
                for my $j (reverse 0 .. 2){
                    if (substr($tqcod, $j, 1) ne "-"){    
                        $trpos -- ;
                    } else {
                        last;
                    }
                }
            } else {
                for my $j (reverse 0 .. 2){
                    if (substr($trcod, $j, 1) ne "-"){
                        $trpos --;
                    } else {
                        last;
                    }
                }
            }
            print "\tENDFS$trpos";
            $is_fs =  0;
            $end_fs = 1;
        }
        my @rbases = split(//, $rcod);
        my @qbases = split(//, $qcod);
        my @out;
        my $basecount = $rpos - 3;
        for my $j (0 .. $#rbases){
            $basecount++;
            if ($rbases[$j] ne $qbases[$j]){
                push @out, "$rbases[$j]$basecount$qbases[$j]";
            }
        }
        if (@out){
            print "\t", join(",", @out) unless $end_fs;
            my ($raa, $qaa) = ($trans{$rcod}, $trans{$qcod});
            if ($raa and $qaa){
                if ($raa ne $qaa){
                    my $apos = $rpos/3;
                    print ",[$raa$apos$qaa]";
                }
            }
        }
        ($rcod, $qcod) = ("") x 2;
        next;
    } elsif (length($rcod) == 3 or length($qcod) == 3) {
        $rcod = "" if length($rcod) == 3;
        $qcod = "" if length($qcod) == 3;
        next if ($is_fs);
        $is_fs = 1;
        my $trcod = substr($refseq, length($refseq)-3, 3);
        my $tqcod = substr($qseq, length($refseq)-3, 3);
        my $trpos = $rpos - 3;
        if ($tqcod =~ m/-/){
            for my $j (0 .. 2){
                if (substr($tqcod, $j, 1) eq "-"){
                    $trpos += $j;
                    last;
                }
            }
        } else {
            for my $j (0 .. 2){
                if (substr($trcod, $j, 1) eq "-"){
                    $trpos += $j;
                    last;
                }
            }
        }
        print "\tFS$trpos";
    }
    print "\n";
    
}


sub translate {
    #table is gencode 11 taken from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
    my $aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    my $b1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
    my $b2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
    my $b3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
    my @aa_a = split("", $aa);
    my @b1_a = split("", $b1);
    my @b2_a = split("", $b2);
    my @b3_a = split("", $b3);
    my %out;
    for my $i (0 .. $#aa_a){
        $out{"$b1_a[$i]$b2_a[$i]$b3_a[$i]"} = $aa_a[$i];
    }
    return (%out);
}

sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}
