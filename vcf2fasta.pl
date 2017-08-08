use List::MoreUtils qw(any all uniq);
use POSIX qw(ceil floor);
my $gff;
my $gzip=0;
my $ref;
my $vcf;
my $phased=0;
my $wref=0;
my %REF=();
my %DATA=();
my %GENE=();
my @lines=();
my @ind;
my $feat="CDS";
my $usage = "Usage:
perl vcf2fasta.pl -f <fasta-ref> -v <vcf-file> -g <gff-file> -e <gff-feature> [ --ref ] [ --phased ]

option --ref will include the reference sequence.
option --phase informs the program that the VCF is pased.
defaults ignore the two options.

examples:
perl vcf2fasta.pl -f ref.fas -v snps.vcf -g annotation.gff -e CDS
perl vcf2fasta.pl -f ref.fas -v snps.phased.vcf -g annotation.gff -e CDS --phased
perl vcf2fasta.pl -f ref.fas -v snps.phased.vcf -g annotation.gff -e CDS --phased --ref
";

if ((scalar(@ARGV) == 0) or (grep { /-+he{0,1}l{0,1}p{0,1}/ } @ARGV)){
     die $usage;
}

if (my ($indV) = grep { $ARGV[$_] =~ /^-v$/ } 0 .. $#ARGV){
     $vcf = $ARGV[$indV+1];
     if ($vcf =~ m/\.gz$/){
          $gzip = 1;
     }
} else {
     die $usage;
}

if (my ($indE) = grep { $ARGV[$_] =~ /^-e$/ } 0 .. $#ARGV){
     $feat = $ARGV[$indE+1];
} else {
     die $usage;
}

if (my ($indG) = grep { $ARGV[$_] =~ /^-g$/ } 0 .. $#ARGV){
     $gff = $ARGV[$indG+1];
} else {
     die $usage;
}

if (my ($indP) = grep { $ARGV[$_] =~ /^--phased$/ } 0 .. $#ARGV){
     $phased = 1;
}

if (my ($indR) = grep { $ARGV[$_] =~ /^--ref$/ } 0 .. $#ARGV){
     $wref= 1;
}

if (my ($indF) = grep { $ARGV[$_] =~ /^-f$/ } 0 .. $#ARGV){
     $ref = $ARGV[$indF+1];
     open(F,"<",$ref) or die "Cannot open $ref.\n";
     print "\n Reading FASTA reference...";
     while(<F>){
          next if (/^\s*$/);
          chomp($_);
          $_ =~ s/\s//g;
          push @lines, $_;
     }
     print " done\n";
} else {
     die $usage;
}

my @indH = grep { $lines[$_] =~ /^>/ } 0 .. $#lines;

print " " . scalar(@indH) . " scaffolds/chromosomes found\n";

print "\n **************\n\n";

for my $i (0 .. $#indH){
     if ($i != $#indH){
          $REF{substr($lines[$indH[$i]],1)} = join('', @lines[ $indH[$i]+1 .. $indH[$i+1]-1 ]);
          $REF{substr($lines[$indH[$i]],1)} =~ tr/a-z/A-Z/;
     } else {
          $REF{substr($lines[$indH[$i]],1)} = join('', @lines[ $indH[$i]+1 .. $#lines ]);
          $REF{substr($lines[$indH[$i]],1)} =~ tr/a-z/A-Z/;
     }
}

print " Converting mVCF to sequence data...\n";

my $count=0;
if ($gzip == 1){
     open(V,"-|","gunzip < $vcf") or die "Cannot open $vcf\n";
} else {
     open(V,"<",$vcf) or die "Cannot open $vcf\n";
}
local $| = 1;
my $T1 = time;
while(<V>){
     if (/^##/){
          next;
     }
     elsif (/^#CHROM/){
          chomp($_);
          @HEAD = split /\t/, $_;
          ($format_ind) = grep { @HEAD[$_] =~ m/FORMAT/ } 0 .. $#HEAD;
          @ind = @HEAD[ ($format_ind+1) .. $#HEAD ];
          if ($phased == 1){
          	my @tmp = ();
          	map { push @tmp, $_ . "_a"; push @tmp, $_ . "_b" } @ind;
          	@ind = @tmp;
          	print " data is phased\n";
          }
          print " " . (scalar(@ind));
          if ($phased == 1){
          	print " haplotypes found\n";
          } else {
          	print " individuals found\n";
          }
          for $SCAFF (keys %REF){
               if ($wref == 1){
                    $DATA{$SCAFF}{'REF'} = $REF{$SCAFF};
               }
               for $IND (@ind){
                    $DATA{$SCAFF}{$IND} = $REF{$SCAFF};
               }
          }
          print " Done replicating reference to indiv. sequences:\n";
          print " @ind\n";
          print " Converting alleles (this might take a while):\n";
          if ($wref == 1){
          	unshift @ind, 'REF';
          }
     } else {
          chomp($_);
          my @field = split /\t/, $_;
          my @i = ($format_ind+1) .. $#HEAD;
          map { s/:.*// } @field[ @i ];
          foreach(@i){
          	$i = $_;
          	if ($phased == 1){
          		my @ii = split /\|/, $field[$i];
          		my $allele1 = getallele(\$field[3], \$field[4], \$ii[0]);
          		my $allele2 = getallele(\$field[3], \$field[4], \$ii[1]);
          		if ($allele1 eq "*" or $allele2 eq "*" ){
          			die "\n failed al pos: $field[0], $field[1], $allele1|$allele2\n";
          		}
          		eval { substr($DATA{$field[0]}{"$HEAD[$i]_a"}, ($field[1]-1), length($allele1), $allele1) };
          		eval { substr($DATA{$field[0]}{"$HEAD[$i]_b"}, ($field[1]-1), length($allele2), $allele2) };
          		die "\n failed al pos: $field[0], $field[1], $allele\n" if $@;
          	} else {
          		if ($field[$i] =~ m/\|/){
          			die "Must specify the --phased flag for phased data\n";
          		}
          		my $allele = getallele(\$field[3], \$field[4], \$field[$i]);
          		if ($allele eq "N"){
          			die "\n failed al pos: $field[0], $field[1], $allele\n";
          		}
          		eval { substr($DATA{$field[0]}{$HEAD[$i]}, ($field[1]-1), length($allele), $allele) };
          		die "\n failed al pos: $field[0], $field[1], $allele\n" if $@;
          	}
          }
          $c = ++$count;
          printf(" SNP site # %10s, chr: %15s, pos: %10s\r", $c, $field[0], $field[1]);
     }
}
local $| = 0;
print "\n All alleles converted\n";
my $T2 = time;
my $T3 = $T2-$T1;
timer($T3);

print "\n **************\n\n";

print " Looking for features...\n";

my $count=0;
my $seq;
local $| = 1;
open(G,"<",$gff) or die "Cannot open $gff\n";
while(<G>){
     next if (/^#/);
     chomp($_);
     my @field = split(/\t/, $_);
     if ($DATA{$field[0]}){
          if ($field[2] eq $feat){
               my ($geneid) = $field[8] =~ m/[= ]"{0,1}(.+?)"{0,1} {0,1};/;
               if (length($geneid) == 0){
               	  $geneid = $field[8];
               }
               foreach(@ind){
                    my $start = $field[3]-1;
                    my $len = ($field[4]-$field[3])+1;
                    $seq = substr($DATA{$field[0]}{$_},$start,$len);
                    $GENE{$_}{$geneid}{'seq'} .= $seq;
                    $GENE{$_}{$geneid}{'sen'} = $field[6];
               }
          }
          $c = ++$count;
          print " GFF line: $c, chr: $field[0], pos: $field[3]               \r";
     }
}
local $| = 0;
my $T4 = time;
my $T5 = $T4-$T2;
print "\n" and timer($T5);

print "\n **************\n\n";

unless (-d $feat){
     `mkdir $feat`;
}

print " Printing features to file...\n";

for $gene (sort {$a cmp $b} keys %{$GENE{$ind[0]}}){
     open(OUT,">","$feat/$gene");
     for $ind (@ind){
          $seq = $GENE{$ind}{$gene}{'seq'};
          if ($GENE{$ind}{$gene}{'sen'} eq '-'){
               $seq = revcomp($seq);
          }
          print OUT ">$ind\n$seq\n";
     }
}
my $T6 = time;
my $TT = $T6-$T1;
print " Total:\n";
timer($TT);

sub revcomp {
     my ($dna) = @_;
     my $rcdna = reverse($dna);
     $rcdna =~ tr/ACTGactg/TGACtgct/;
     return($rcdna);
}

sub getallele {
     my ($ref,$alt,$gen) = @_;
     my %trans=();
     $trans{0} = $$ref;
     $trans{"."} = "?";
     if ($$alt =~ m/.,./){
          my @al = split(/,/, $$alt);
          for my $i (0 .. $#al){
               $trans{$i+1} = $al[$i];
          }
     } else {
          $trans{1} = $$alt;
     }
     if ($$gen =~ m/.\/./){
          my @both = split /\//, $$gen; 
          my $Ap = getiupac(map { $_ = $trans{$_} } @both);
          return($Ap);
     }
     elsif ($$gen =~ m/.{1}/){
          return($trans{$$gen});
     } else {
     	return('N');
     }
}

sub timer {
     my ($T) = @_;
     if ($T < 86400){
          if ($T > 3600){
               $T = sprintf("%.2f",$T/3600);
               print " Done in $T hours\n";
          }
          elsif ($T > 60){
               $T = sprintf("%.2f",$T/60);
               print " Done in $T minutes\n";
          }
          else {
               print " Done in $T seconds\n";
          } 
     }
     else {
          $T = sprintf("%.2f",$T/86400);
          print " Done in $T days\n";
     }
}

sub getiupac {
     my @als = @_;
     if (scalar(grep { /A/ } @als) == 2){
          return 'A';
     }
     elsif (scalar(grep { /T/ } @als) == 2){
          return 'T';
     }
     elsif (scalar(grep { /C/ } @als) == 2){
          return 'C';
     }
     elsif (scalar(grep { /G/ } @als) == 2){
          return 'G';
     }
     elsif (scalar(grep { /A/ } @als) == 1 and scalar(grep { /T/ } @als) == 1){
          return 'W';
     }
     elsif (scalar(grep { /A/ } @als) == 1 and scalar(grep { /C/ } @als) == 1){
          return 'M';
     }
     elsif (scalar(grep { /A/ } @als) == 1 and scalar(grep { /G/ } @als) == 1){
          return 'R';
     }
     elsif (scalar(grep { /T/ } @als) == 1 and scalar(grep { /C/ } @als) == 1){
          return 'Y';
     }
     elsif (scalar(grep { /T/ } @als) == 1 and scalar(grep { /G/ } @als) == 1){
          return 'K';
     }
     elsif (scalar(grep { /G/ } @als) == 1 and scalar(grep { /C/ } @als) == 1){
          return 'S';
     }
     else{
          return '?';
     }
}

