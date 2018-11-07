#!/usr/bin/perl
use warnings;
use strict; 
use Bio::DB::GenBank;
use Bio::SeqIO;  
use Text::Wrap;
use Data::Dumper;
#use lib '/work/home/OM/wangj/project/perl'; 
#use BioSeq;

die "\nUsage: perl neoantigen.pl HLA_file mutation_file output_file sample_id_col gene_id_col NM_id_col cdot_col pdot_col vaf_col aa_length\n\n" unless @ARGV==10;

`mkdir -p Peptide`;
`mkdir -p NetMHC`;

#my $wp =  get_fsCDS("NM_006572","c.788dup","p.N263Kfs*11");
#print "$wp\n";

open HLA,$ARGV[0] or die $!;
open CP,$ARGV[1] or die $!;
open NA,">$ARGV[2]" or die $!;


my $id_col = int($ARGV[3]) - 1;
my $gene_col = int($ARGV[4]) - 1;
my $nm_col = int($ARGV[5]) - 1;
my $cdot_col = int($ARGV[6]) - 1;
my $pdot_col = int($ARGV[7]) - 1;
my $vaf_col = int($ARGV[8]) - 1;
my $aalen = int($ARGV[9]);
my $llen = $aalen + 1;


my %hashHLA;

while (<HLA>){
	chomp;
	my @line =  split(/\t/,$_);
	#my $id = substr($line[0],0,10);
	my $id = $line[0];
	#print "@line[1..6]\n";HLA-C12:02 <- A*31:01
	#map {$_=~ s/([ABC])\*/HLA-$1/} @line[1..6];
	#print "@hla\n";
	@{$hashHLA{$id}}=unique(@line[1..2]);
}
close HLA;

#print Dumper(\%hashHLA);

print NA "ORDER_ID\tGENE\tNM_ID\tCdot\tPdot\tVAF\tAlleles\tAllele\tMT_PEP\tAffinity(nm)\t%Rank\tWT_PEP\tAffinity(nm)\t%Rank\n";

while (<CP>){
	chomp;
	next if /^ORDER_ID|Sample_ID|OM_ID/;
	my @line =  split(/\t/,$_);
	my $id = $line[$id_col];
	my $gene = $line[$gene_col];
	my $nm = $line[$nm_col];
	my $cdot = $line[$cdot_col];
	my $pdot = $line[$pdot_col];
	my $vaf = $line[$vaf_col];
	
	next if ($pdot =~ /[pP]\.[A-Z]\d+\*$/);
	next if ($pdot =~ /\?/);
	next if ($pdot =~ /[pP]\.([A-Z])\d+\1$/);
	next if ($pdot eq "");
	next if ($pdot eq "-");
	next if ($pdot eq "NA");

	my ($mt_pep,$wt_pep);
	
	eval { get_pep($nm,$cdot,$pdot,$aalen);}; # get error

	if ($@){
                $mt_pep = undef;
                $wt_pep = undef;
                print STDERR "$id\t$gene\t$nm\t$cdot\t$pdot have an error. Error: $@\n";
        }else{
                ($mt_pep,$wt_pep) = get_pep($nm,$cdot,$pdot,$aalen);
        }  
	
	my ($mt_out,$wt_out);
	if (defined $mt_pep && defined $wt_pep)	{
		$mt_out = creat_fasta($nm,$pdot,$mt_pep,"./Peptide/${nm}_${pdot}_MT.fasta");
		$wt_out = creat_fasta($nm,$pdot,$wt_pep,"./Peptide/${nm}_${pdot}_WT.fasta");
		
	}elsif(defined $mt_pep && !defined $wt_pep){
		$mt_out = creat_fasta($nm,$pdot,$mt_pep,"./Peptide/${nm}_${pdot}_MT.fasta");
		$wt_out = undef;
	}else{
		$mt_out = undef;
		$wt_out = undef;
	}

	if (defined $mt_out && defined $wt_out){
		for my $alle (@{$hashHLA{$id}}){
			my $mt_mhc_out = run_netmhcII($mt_out,$alle,"./NetMHC/${nm}_${pdot}_MT_${alle}");
			my $wt_mhc_out = run_netmhcII($wt_out,$alle,"./NetMHC/${nm}_${pdot}_WT_${alle}");
			
			#print "$mt_mhc_out\n$wt_mhc_out\n";

			open MO,"$mt_mhc_out" or die "cannot open the file: $mt_mhc_out. $!\n";
			open WO,"$wt_mhc_out" or die "cannot open the file: $wt_mhc_out. $!\n";

			my @ML = <MO>;
			my @WL = <WO>;

			close MO;
			close WO;

			for (my $i=0;$i<=$#ML;$i++){
				#print "$i\n";
				next if $ML[$i] =~ /^#/;
				if ($ML[$i] =~ /\s+DRB1_\d+/ && $ML[$i] =~ /\s+[A-Z]{$llen}/){
					#print "ok\n";
					my @mline = split(/\s+/,$ML[$i]);
					my @wline = split(/\s+/,$WL[$i]);

					print NA "$id\t$gene\t$nm\t$cdot\t$pdot\t$vaf\t",join(";",@{$hashHLA{$id}}),"\t$alle\t$mline[3]\t$mline[9]\t$mline[10]\t$wline[3]\t$wline[9]\t$wline[10]\n";
				}	
			}
		}
	}elsif(defined $mt_out && !defined $wt_out){
		for my $alle (@{$hashHLA{$id}}){
			my $mt_mhc_out = run_netmhcII($mt_out,$alle,"./NetMHC/${nm}_${pdot}_MT_${alle}");
			open MO,"<",$mt_mhc_out or die "cannot open the file: $mt_mhc_out. $!\n";
			my @ML = <MO>;
			close MO;
			for (my $i=0;$i<=$#ML;$i++){
				next if $ML[$i] =~ /^#/;
				if ($ML[$i] =~ /\s+DRB1_\d+/ && $ML[$i] =~ /\s+[A-Z]{$llen}/){
					my @mline = split(/\s+/,$ML[$i]);
					print NA "$id\t$gene\t$nm\t$cdot\t$pdot\t$vaf\t",join(";",@{$hashHLA{$id}}),"\t$alle\t$mline[3]\t$mline[9]\t$mline[10]\tNAN\tNAN\tNAN\n";
				}
			}
		}
	}else{
		for my $alle (@{$hashHLA{$id}}){
			print NA "$id\t$gene\t$nm\t$cdot\t$pdot\t$vaf\t",join(";",@{$hashHLA{$id}}),"\t$alle\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\n";
		}
	}
		

}
close CP;
close NA;


## sub 

sub unique {
	my @array = @_;
	my (%hash, @unique);
	foreach my $a (@array) {
		if (!defined($hash{$a})){
			push(@unique,$a);
			$hash{$a} = 1;
		}
	}
	return @unique ;
}


sub creat_fasta{

	my ($acc,$pdot,$pseq,$outfile) = @_;
	
	open OUT,">$outfile" or die $!;

	print OUT ">${acc}_${pdot}\n$pseq\n";

	return $outfile;
		
}


sub get_CDS { # get peptide from substitution
	my ($acc,$pos,$len,$mut,$mtout,$wtout)=@_;

	open MT,">",$mtout or die "# cannot read $mtout:$!\n";
	open WT,">",$wtout or die "# cannot read $wtout:$!\n";
	#print "$acc,$pos,$len\n";
	
	my $gb   = new Bio::DB::GenBank;
	my $seq1 = $gb->get_Seq_by_acc($acc);
	#my $sequence = $seq1->seq;
	#print "$sequence\n";
	my ($wild_pep,$mut_pep,$aa);

	for my $feat ($seq1->get_SeqFeatures){
		if ($feat->primary_tag eq 'CDS'){
			#print Dumper($feat);
			#print  $feat->get_tag_values('translation'),"\n";
			my @pep = $feat->get_tag_values("translation");
			#print Dumper(\@pep);
			my $prot = $pep[0];
			my ($left,$right);
			if ($pos <= $len ){
				$aa = substr($prot,$pos-1,1);
				$left = substr($prot,0,$pos-1);
				$right = substr($prot,$pos,$len);
			}else{
				$aa = substr($prot,$pos-1,1);
				$left = substr($prot,$pos-$len-1,$len);
				$right = substr($prot,$pos,$len);
			}
			$wild_pep = $left . $aa . $right;
			$mut_pep = $left . $mut . $right;
		}
	}
	
	print MT ">${acc}_p\.${aa}$pos${mut}_MT\n$mut_pep\n";
	print WT ">${acc}_p\.${aa}$pos${mut}_WT\n$wild_pep\n";
	
	return $mtout,$wtout;
}

sub get_fsCDS { # get peptide from frameshift
	my ($acc,$cdot,$pdot)=@_;
	my ($cds,$new_cds,$wt_pep);
	
	my $gb   = new Bio::DB::GenBank;
	my $seq1 = $gb->get_Seq_by_acc($acc);
	my $sequence = $seq1->seq;
	
	for my $feat ($seq1->get_SeqFeatures){
		if ($feat->primary_tag eq 'CDS'){
		#print Dumper($feat);
		my $start = $feat->start;
		my $len   = $feat->length;
		#print "$len\n";
		$cds = substr($sequence,$start-1,$len);
		#print "$cds\n";
		}
	}

	#print "$cds\n";
	if ($cdot =~ /[cC]\.(\d+)dup/){ # dup of one base
		my $base = substr($cds,$1-1,1);
		my $left = substr($cds,0,$1);
		my $right = substr($cds,$1);
		$new_cds = $left . $base . $right;
	}
	#print "$new_cds\n";


	my $seq_obj = Bio::Seq->new(
		-seq => $new_cds,
		-alphabet => 'dna' );	

	my $prot_obj = $seq_obj->translate(-orf => 1, -complete => 1); 
	my $pep_seq =  $prot_obj->seq;

	if ($pdot =~ /[pP]\.([A-Z]+)(\d+)([A-Z]+)fs\*(\d+)/){
		print "$2\n";
		$wt_pep = substr($pep_seq,$2-9);
	}
	
	return $wt_pep;
}

sub run_netmhc {
	
	my ($pep,$len,$allele,$outfile) = @_ ;

	#my $netMHC = '/work/home/OM/wangj/program/HLA/netMHCIIpan-3.2/netMHCIIpan';
	my $netMHC = '/work/home/OM/wangj/program/HLA/netMHCpan-4.0/netMHCpan';

        system("$netMHC -l $len -a $allele -BA '$pep' > '$outfile' ");

	return $outfile;
}

sub run_netmhcII {
	
	my ($pep,$allele,$outfile) = @_ ;
	
	my $netMHCII = '/work/home/OM/wangj/program/HLA/netMHCIIpan-3.2/netMHCIIpan';

	system("$netMHCII -f '$pep' -a $allele > '$outfile'");
	
	return $outfile;

}

sub get_pep {
	
	my ($acc,$cdot,$pdot,$len)=@_;

	my $gb   = new Bio::DB::GenBank;
        my $seq1 = $gb->get_Seq_by_acc($acc);
        my $sequence = $seq1->seq;
        #print "$sequence\n";
		
	my ($prot,$cds,$mrna);
	for my $feat ($seq1->get_SeqFeatures){
                if ($feat->primary_tag eq 'CDS'){
                        #print Dumper($feat);
                        #print  $feat->get_tag_values('translation'),"\n";
                        my @pep = $feat->get_tag_values("translation");
                        #print Dumper(\@pep);
                        $prot = $pep[0];

			my $start = $feat->start;
	                my $len   = $feat->length;
        	        #print "$len\n";
			#$cds = substr($sequence,$start-1,$len);
			$mrna = substr($sequence,$start-1);
			$cds = $mrna;
                }
        }
	
        
	my ($mtpep,$wtpep);

	if ($pdot =~ /[pP]\.([A-Z])(\d+)([A-Z])$/){# substitution
		my $wp = $1; my $pos = $2; my $mp = $3;
		my $aa = substr($prot,$pos-1,1);
		my ($left,$right);
		if ($aa ne $wp ){
			print STDERR "Error: $acc $pdot have an error transcript.\n";
			return undef,undef;
		}else{
			if ($pos <= $len){
				$left = substr($prot,0,$pos-1);
                                $right = substr($prot,$pos,$len);
			}else{
				$left = substr($prot,$pos-$len-1,$len);
				$right = substr($prot,$pos,$len);
			}	
		}
		
		$mtpep = $left . $mp . $right;
		$wtpep = $left . $aa . $right;	
	}elsif($pdot =~ /[pP]\.([A-Z])(\d+)([A-Z])fs\*\d+/){ # frame shift and translation termination codon (stop codon, no-stop change)
		my $wp = $1; my $pos = $2; my $mp = $3;
		my $new_cds;
		if ($cdot =~ /[cC]\.(\d+)[ATCG]\>([ATCG])/){ # substitution on the stop codon
			substr($cds,$1-1,1) = $2;
			$new_cds = $cds;
		
		}elsif($cdot =~ /[cC]\.(\d+)dup[A-Z]?/) { # cdot dup of one base
			my $base = substr($cds,$1-1,1);
                	my $left = substr($cds,0,$1);
	                my $right = substr($cds,$1);
        	        $new_cds = $left . $base . $right;
		}elsif($cdot =~ /[cC]\.(\d+)_(\d+)dup[A-Z]*/){ # cdot dup of more than one base
			my $base = substr($cds,$1-1,$2-$1+1);
			my $left = substr($cds,0,$2);
			my $right = substr($cds,$2);
			$new_cds = $left . $base . $right;
		}elsif($cdot =~ /[cC]\.(\d+)del[A-Z]?/){ # cdot del of one base
			substr($cds,$1-1,1) = "";
			$new_cds = $cds;			
		}elsif($cdot =~ /[cC]\.(\d+)_(\d+)del[A-Z]*/){ # cdot del of more than one base
			substr($cds,$1-1,$2-$1+1) = "";
			$new_cds = $cds;			
		}elsif($cdot =~ /[cC]\.(\d+)_(\d+)ins([A-Z]+)/){ # cdot ins 
			$new_cds = substr($cds,0,$1) . $3 . substr($cds,$2-1);
		}elsif($cdot =~ /[cC]\.(\d+)_(\d+)delins([A-Z]+)/){
			substr($cds,$1-1,$2-$1+1) = $3;
			$new_cds = $cds;			
		}else{
			return undef,undef;
		}
		
		my $seq_obj = Bio::Seq->new(
 	               -seq => $new_cds,
        	        -alphabet => 'dna' );

        	my $prot_obj = $seq_obj->translate(-orf => 1, -complete => 1);
	        my $pep_seq =  $prot_obj->seq;

		if ($pos <= $len){
			$mtpep =  $pep_seq ;
		}else{
			$mtpep = substr($pep_seq,$pos-$len-1); 
		}

		$wtpep = undef;
	}elsif($pdot =~ /[pP]\.(?:\*|Ter)(\d+)([A-Z])ext\*\d+/){ # translation termination codon (stop codon, no-stop change) 
                #my $wp = $1; 
                my $pos = $1; my $mp = $2;
                my $new_cds;
                if ($cdot =~ /[cC]\.(\d+)[ATCG]\>([ATCG])/){ # substitution on the stop codon
                        #print "$cds\n";
                        substr($mrna,$1-1,1) = $2;
                        $new_cds = $mrna;
                        #print "$new_cds\n";

                }elsif($cdot =~ /[cC]\.(\d+)dup[A-Z]?/) { # cdot dup of one base
                        my $base = substr($mrna,$1-1,1);
                        my $left = substr($mrna,0,$1);
                        my $right = substr($mrna,$1);
                        $new_cds = $left . $base . $right;
                }elsif($cdot =~ /[cC]\.(\d+)_(\d+)dup[A-Z]*/){ # cdot dup of more than one base
                        my $base = substr($mrna,$1-1,$2-$1+1);
                        my $left = substr($mrna,0,$2);
                        my $right = substr($mrna,$2);
                        $new_cds = $left . $base . $right;
                }elsif($cdot =~ /[cC]\.(\d+)del[A-Z]?/){ # cdot del of one base
                        substr($mrna,$1-1,1) = ""; 
                        $new_cds = $mrna;    
                }elsif($cdot =~ /[cC]\.(\d+)_(\d+)del[A-Z]*/){ # cdot del of more than one base
                        substr($mrna,$1-1,$2-$1+1) = ""; 
                        $new_cds = $mrna;    
                }elsif($cdot =~ /[cC]\.(\d+)_(\d+)ins([A-Z]+)/){ # cdot ins 
                        $new_cds = substr($mrna,0,$1) . $3 . substr($mrna,$2-1);
                }elsif($cdot =~ /[cC]\.(\d+)_(\d+)delins([A-Z]+)/){
                        substr($mrna,$1-1,$2-$1+1) = $3;
                        $new_cds = $mrna;
                }else{
                        return undef,undef;
                }

                my $seq_obj = Bio::Seq->new(
                       -seq => $new_cds,
                        -alphabet => 'dna' );

                my $prot_obj = $seq_obj->translate(-orf => 1, -complete => 1);
                my $pep_seq =  $prot_obj->seq;

                if ($pos <= $len){
                        $mtpep =  $pep_seq ;
                }else{
                        $mtpep = substr($pep_seq,$pos-$len-1);
                }

                $wtpep = undef;
	
	}elsif($pdot =~ /[pP]\.([A-Z])(\d+)del/){ # del of one base
		my $wt = $1; my $pos = $2;
		if ($pos <= $len){
			$mtpep = substr($prot,0,$pos-1) . substr($prot,$pos,$len);
		}else{
			$mtpep = substr($prot,$pos-$len-1,$len) . substr($prot,$pos,$len);
		}
		
		$wtpep = undef;
	}elsif($pdot =~ /[pP]\.[A-Z](\d+)_[A-Z](\d+)del/){ # del of more than one base
		my $pos1 = $1; my $pos2 = $2;
		if ($pos1 <= $len){
			$mtpep = substr($prot,0,$pos1-1) . substr($prot,$pos2,$len);
		}else{
			$mtpep = substr($prot,$pos1-$len-1,$len) . substr($prot,$pos2,$len);	
		}
		
		$wtpep = undef;
	}elsif($pdot =~ /[pP]\.([A-Z])(\d+)dup/){ # dup of one base
		my $wt = $1; my $pos = $2;
		if ($pos < $len){
			$mtpep = substr($prot,0,$pos) . $wt . substr($prot,$pos,$len);
		}else{
			$mtpep = substr($prot,$pos-$len,$len) . $wt . substr($prot,$pos,$len);
		}
		
		$wtpep = undef;
	}elsif($pdot =~ /[pP]\.[A-Z](\d+)_[A-Z](\d+)dup/){ # dup of more than one base
		my $pos1 = $1; my $pos2 = $2;
		my $base = substr($prot,$pos1-1,$pos2-$pos1+1);
		if ($pos2 < $len){
			$mtpep = substr($prot,0,$pos2) . $base . substr($prot,$pos2,$len); 	
		}else{
			$mtpep = substr($prot,$pos2-$len,$len) . $base . substr($prot,$pos2,$len);
		}

		$wtpep = undef;
	}elsif($pdot =~ /[pP]\.[A-Z](\d+)_[A-Z](\d+)ins([A-Z]+)/){ # ins
		my $pos1 = $1; my $pos2 = $2; my $base = $3;
		if ($pos1 < $len){
			$mtpep = substr($prot,0,$pos1) . $base . substr($prot,$pos2-1,$len);
		}else{
			$mtpep = substr($prot,$pos1-$len,$len) . $base . substr($prot,$pos2-1,$len);
		}
		
		$wtpep = undef;
	}elsif($pdot =~ /[pP]\.[A-Z](\d+)_[A-Z](\d+)delins([A-Z]+)/){ #delins
		my $pos1 = $1; my $pos2 = $2; my $base = $3;
		if ($pos1 <= $len){
			$mtpep = substr($prot,0,$pos1-1) . $base . substr($prot,$pos2,$len);
		}else{
			$mtpep = substr($prot,$pos1-$len-1,$len) . $base . substr($prot,$pos2,$len);
		}

		$wtpep = undef;	
	}else{
		$mtpep = undef;
		$wtpep = undef;
	}

	return $mtpep,$wtpep;

}





