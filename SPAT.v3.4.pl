#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Array::Utils qw(:all);
use List::Util qw(max);
use List::MoreUtils qw(firstidx);


###  Reading options

my ($help, $man, $input, $format, $best_score, $graph, $score_threshold, $output_path, $update);

GetOptions(
	'help|?' => \$help,
	man => \$man,
	"input=s" => \$input,
	"format=s" => \$format,
	"best_score=s" => \$best_score,
	"graph=s" => \$graph,
	"score_threshold=i" => \$score_threshold,
	"output=s" => \$output_path,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


=head1 NAME
 
SPAT (Surface Protein Annotation Tool) Version 3.4 - updated jan-2023

=head1 AUTHOR

Jean-Francois Spinella, E<lt>jfspinella@gmail.comE<gt>
Team Sauvageau (Institute for Research in Immunology and Cancer)

=head1 SYNOPSIS

>perl SPAT.v3.4.pl -help [brief help message] -man [full documentation]

>annotation: perl SPAT.v3.4.pl -input path/to/your/input/file.txt -format uniprot_id|ensembl_id|gene_name -output [path/to/your/output/directory/]

=head1 DESCRIPTION

>SPAT takes either uniprot_id, ensembl_id or gene_name as input in a .txt file and returns 1 SPAT score (and associated annotations)

>SPAT requires an input file (.txt format) containing either a list of gene_names or a list of uniprot_ids (1 per line)

>Please, don't modify files located in the SPAT/database directory

=head1 DATE 

Jan-2023

=head1 REQUIREMENTS

>The following perl librairies must be installed: 

    -Array::Utils (installation: cpan Array::Utils)

    -Pod::Usage (installation: cpan Pod::Usage)

=head1 OPTIONS

-help <brief help message>

-man <full documentation>

-input <path/to/your/input/file.txt> mandatory

-format <uniprot_id|ensembl_id|gene_name> mandatory

-best_score <true|false> optional

-graph <true|false> optional

-score_threshold <0 to 10> optional

-output <path/to/your/output/directory/> optional     
      
=cut  



### Test for options

print "\n##############################################\n";

if(!defined $input){
    print "\nWARNING: THE RAW DATA FILE (AND ITS COMPLETE PATH) HAS TO BE PROVIDED AS ARGUMENT\n";
    exit 0;
}

if($format eq "uniprot_id" || $format eq "ensembl_id" || $format eq "gene_name"){
}else{
    print "\nWARNING: THE FORMAT OF INPUT HAS TO BE PROVIDED AS uniprot_id, ensembl_id OR gene_name\n";
    exit 0;
}

my $file_name;
if(!defined $output_path){
    if($input =~ /(.*\/)([^\/]+).txt/){
        $output_path = $1;
        $file_name = $2;
    }
    print "\nWARNING: No output path has been provided, the output will be created in the input directory.\n";
}else{
    if($output_path !~ /\/$/){
        $output_path = $output_path."/";
    }
    if($input =~ /(.*\/)([^\/]+).txt/){
        $file_name = $2;
    }
    print "\n- The output will be created in: ".$output_path."\n";
}

if(!defined $best_score){
    print "\n- By default, if several references exist for a protein, only the entry with the best SPAT score will be reporter. Use the option \"-best_score false\" to output all entries.\n";
    $best_score = "true";
}elsif($best_score eq "false"){
    print "\n- The option \"-best_score\" has been turned off. All entries will be reported (several for a protein if several references exist). The SPAT score may differ between entries according to associated annotations.\n";
}else{
    print "\n- The option \"-best_score\" has been turned on. If several references exist for a protein, only the entry with the best SPAT score will be reporter. Use the option \"-best_score false\" to output all entries (default value is true).\n";
    $best_score = "true";
}

if(!defined $score_threshold){
    print "\n- No score threshold has been set (default threshold is >=0).\n";
    $score_threshold = 0;
}elsif($score_threshold >= 0 and $score_threshold <=10){
    print "\n- Score threshold has been set to $score_threshold. Only entries with a SPAT score >= $score_threshold will be reported.\n";
}else{
    print "\n- If used, the option \"-score_threshold\" must take a value between 0 and 10 (default threshold is >=1).\n";
    $score_threshold = 0;
}
    
if(!defined $graph or $graph eq "false"){
    print "\n- Graph option has not been selected or has been turned off, no score distribution graph will be created. The graph can still be created by running the command line provided in the log file.\n\n";
    $graph = "false";
}elsif($graph eq "true"){
    print "\n- The score distribution graph will be created in the same directory.\n\n";
}else{
    print "\n- Use the option \"-graph true\" to create the score distribution graph (default value is false).\n\n";
    $graph = "false";
}


### Define var

my (@S1, @S2, @S3, @S4, @S5, @RF, @SF);
my $datestring = localtime();



### Score db

my $file = "./SPAT_scores.v4.txt";
my %annotation;
open (FILE, "< $file") or die "Can't open $file: $!";
while (my $l =<FILE>){
    chomp $l;
    next if $l =~ /^#/;
    my @l = split(/\t/, $l);
    if($l[2] eq "S1"){
        push @S1, $l[1];
    }elsif($l[2] eq "S2"){
        push @S2, $l[1];
    }elsif($l[2] eq "S3"){
        push @S3, $l[1];
    }elsif($l[2] eq "S4"){
        push @S4, $l[1];
    }elsif($l[2] eq "S5"){
        push @S5, $l[1];
    }elsif($l[2] eq "RF"){
        push @RF, $l[1];
    }elsif($l[2] eq "SF"){
        push @SF, $l[1];
    }
    $annotation{$l[1]}=$l[0];
}

close FILE;


### Pubmed-term db 

my $file2 = "./database/ndb_all.ref.txt";
my %ref;
open (FILE2, "< $file2") or die "Can't open $file2: $!";
while (my $l =<FILE2>){
    chomp $l;
    my @l = split(/\t/, $l);
    $ref{$l[0]."\t".$l[1]}=$l[2];
}

close FILE2;



### Create hash with ndb_all.term info

my $file3 = "./database/ndb_all.term.txt";
my %spat;
open (FILE3, "< $file3") or die "Can't open $file3: $!";
while (my $l =<FILE3>){
    chomp $l;
    my @l = split(/\t/, $l);    
    my $membrane_score = 0;
    my $membrane_class = "NA";
    my @membrane_terms = ();
    my @annotation = ();
    my @ref = ();
    my @rf_terms = ();
    my @sf_terms = ();
    my $S1_count = 0;
    my $S2_count = 0;
    my $S3_count = 0;
    my $S4_count = 0;
    my $S5_count = 0;
    my $RF_count = 0;
    my $SF_count = 0;

    if (exists $l[1]){
            
            my @term = split(/,/, $l[1]);
            $S1_count = scalar(intersect(@S1,@term));
            $S2_count = scalar(intersect(@S2,@term));
            $S3_count = scalar(intersect(@S3,@term));
            $S4_count = scalar(intersect(@S4,@term));
            $S5_count = scalar(intersect(@S5,@term));
            $RF_count = scalar(intersect(@RF,@term));
            $SF_count = scalar(intersect(@SF,@term));
        
        
            ## CALCULATE A MEMBRANE SCORE ##
            
            if($S1_count>0){
                
                $membrane_score=10;
                $membrane_class="S1";
                $membrane_terms[0]="Cell surface receptor";
                push(@membrane_terms, intersect(@S1,@term));
            
            }elsif($S2_count>0){
                
                $membrane_score=9;
                $membrane_class="S2";
                $membrane_terms[0]="Assoc. with external part of cell membrane";
                push(@membrane_terms, intersect(@S2,@term));
                
            }elsif($S3_count>0){
            
                $membrane_score=8;
                $membrane_class="S3";
                $membrane_terms[0]="Potentially assoc. with external part of cell membrane";
                push(@membrane_terms, intersect(@S3,@term));
                
                if($S4_count>0 and $S5_count==0){
                    
                    $membrane_score-=1;;
                    $membrane_class="S3,S4";
                    $membrane_terms[scalar @membrane_terms]="Assoc. with cell membrane";
                    push(@membrane_terms, intersect(@S4,@term));
                    
                    if($S4_count>=$S3_count+2){
                        
                        $membrane_score-=1;
                        $membrane_class="S3,S4+";
                        
                    }
                    
                }elsif($S5_count>0){
                
                    $membrane_score-=2;
                    $membrane_class="S3,S5";
                    $membrane_terms[scalar @membrane_terms]="Assoc. with membrane";
                    push(@membrane_terms, intersect(@S5,@term));
                    
                    if($S5_count>=$S3_count+2){
                        
                        $membrane_score-=1;
                        $membrane_class="S3,S5+";
                        
                    }
                    
                }
            
            }elsif($S4_count>1){
            
                $membrane_score=5;
                $membrane_class="S4";
                $membrane_terms[0]="Assoc. with cell membrane";
                push(@membrane_terms, intersect(@S4,@term));
                
                if($S5_count>0){
                    
                    $membrane_score-=1;
                    $membrane_class="S4,S5";
                    $membrane_terms[scalar @membrane_terms]="Assoc. with membrane";
                    push(@membrane_terms, intersect(@S5,@term));
                            
                    if($S5_count>=$S4_count+2){
                            
                        $membrane_score-=1;
                        $membrane_class="S4,S5+";
                            
                    }
                        
                }
            
            }elsif($S4_count==1){
                
                $membrane_score=3;
                $membrane_class="S4(1)";
                $membrane_terms[0]="Assoc. with cell membrane";
                push(@membrane_terms, intersect(@S4,@term));
                    
                    if($S5_count>0){
                       
                        $membrane_score-=1;
                        $membrane_class="S4(1),S5";
                        $membrane_terms[scalar @membrane_terms]="Assoc. with membrane";
                        push(@membrane_terms, intersect(@S5,@term));
                    
                        if($S5_count>=$S4_count+2){
                            
                            $membrane_score-=1;
                            $membrane_class="S4(1),S5+";
                        
                        }
                    
                    }
                    
            }elsif($S5_count>=1){
                
                $membrane_score=1;
                $membrane_class="S5";
                $membrane_terms[0]="Assoc. with membrane";
                push(@membrane_terms, intersect(@S5,@term));
                        
            }
            
            
            # RED FLAG SCORE MODIF
            if($membrane_class eq "S1" or $membrane_class eq "S2"){
                if($RF_count>1){
                    $membrane_score-=1;
                    $membrane_class = $membrane_class.",RF+";
                }elsif($RF_count==1){
                    $membrane_score-=1;
                    $membrane_class = $membrane_class.",RF";
                }
            }else{
                if($RF_count>1){
                    $membrane_score-=3;
                    $membrane_class = $membrane_class.",RF+";
                }elsif($RF_count==1){
                    $membrane_score-=1;
                    $membrane_class = $membrane_class.",RF";
                }
            }
            if($RF_count>0){
                push(@rf_terms, intersect(@RF,@term));
            }else{
                @rf_terms = ("NA");
            }
            
            
            # SECRETED FLAG 
            if($membrane_class =~ /S1/ or $membrane_class =~ /S2/){
                if($SF_count>1){
                    $membrane_class = $membrane_class.",SF+";
                }elsif($SF_count==1){
                    $membrane_class = $membrane_class.",SF";
                }
            }else{
                if($SF_count>1){
                    $membrane_score-=2;
                    $membrane_class = $membrane_class.",SF+";
                }elsif($SF_count==1){
                    #$membrane_score-=1;
                    $membrane_class = $membrane_class.",SF";
                }
            }
            if($SF_count>0){
                push(@sf_terms, intersect(@SF,@term));
            }else{
                @sf_terms = ("NA");
            }

                        
            # ANNOTATION
            if($membrane_score == 0){
                @annotation = ("NA");
            }else{
                foreach my $membrane_terms(@membrane_terms){
                    if(exists $annotation{$membrane_terms}){
                        push(@annotation, $annotation{$membrane_terms});
                    }
                }
            }
            
            
            # PROT WITH NO TERM MATCHING THE SURF LIST
            if($S1_count+$S2_count+$S3_count+$S4_count+$S5_count==0){
                $membrane_score=0;
                @membrane_terms = ("NA");
            }else{
                if($membrane_score<1){
                    $membrane_score=1;
                }
            }
            
            
            # PUBMED REF
            if($membrane_score == 0){
                @ref = ("NA");
            }else{
                foreach my $membrane_terms(@membrane_terms){
                    if(exists $ref{$l[0]."\t".$membrane_terms}){
                        push(@ref, $membrane_terms.">".$ref{$l[0]."\t".$membrane_terms});
                    }
                }
            }
            if(scalar @ref == 0){
                @ref = ("NA");
            }
            
        }
        
        
        ## CREATE HASH ##
        $spat{$l[0]} = "$membrane_score\t".join( ',', @ref )."\t$RF_count\t".join( ',', @rf_terms )."\t$SF_count\t".join( ',', @sf_terms )."\t".$l[2];
       
}

close FILE3;



### Annotate raw data and create output

# translator
my $file4 = "./database/translator.txt";
my %transl; my %gene_name; my %ensembl_id;
open (FILE4, "< $file4") or die "Can't open $file4: $!";
while (my $l =<FILE4>){
    chomp $l;
    my @l = split(/\t/, $l);
    if(exists $spat{$l[0]}){
        if($format eq "uniprot_id"){
            my @ndb_id;
            if(exists $transl{$l[1]}){
                # test for redundancy in ndb ids
                my $test = grep(/$l[0]/i, @{$transl{$l[1]}});
                if($test > 0){
                    next;
                }else{
                    push( @{$transl{$l[1]}}, $l[0]); 
                }
            }else{
                $ndb_id[0] = $l[0];
                $transl{$l[1]} = [ @ndb_id ];
            }
        }elsif($format eq "ensembl_id"){
            my @ndb_id;
            if(exists $transl{$l[2]}){
                # test for redundancy in ndb ids
                my $test = grep(/$l[0]/i, @{$transl{$l[2]}});
                if($test > 0){
                    next;
                }else{
                    push( @{$transl{$l[2]}}, $l[0]); 
                }
            }else{
                $ndb_id[0] = $l[0];
                $transl{$l[2]} = [ @ndb_id ];
            }
        }elsif($format eq "gene_name"){
            my @ndb_id;
            if(exists $transl{$l[3]}){
                # test for redundancy in ndb ids
                my $test = grep(/$l[0]/i, @{$transl{$l[3]}});
                if($test > 0){
                    next;
                }else{
                    push( @{$transl{$l[3]}}, $l[0]); 
                }
            }else{
                $ndb_id[0] = $l[0];
                $transl{$l[3]} = [ @ndb_id ];
            }
        }
        $gene_name{$l[0]} = $l[3];
        $ensembl_id{$l[2]} = $l[0];
    }
}
close FILE4;

# human prot atlas
my $file5 = "./database/prot_atl.txt";
my %atlas;
my @prognostics_header;
my @expression_header;

open (FILE5, "< $file5") or die "Can't open $file5: $!";
while (my $l =<FILE5>){
    
    chomp $l;
    $l =~ tr/"//d;
    my @l = split(/\t/, $l);
    my @output = ("NA","NA","NA","NA","NA");
    
    # main location
    if($.>1 and $l[53] ne ""){
        $output[0] = $l[53];
    }
    
    # antibody
    if($.>1 and $l[55] ne ""){
        my @antibody;
        while ($l[55] =~ m/\w+\d+: (AB_\d+)/g) {
            push @antibody,$1;
        }
        $output[1] = join(",",@antibody);
    }
    
    # prot level
    if($.>1 and $l[227] ne ""){
        $output[4] = $l[227];
    }elsif($.>1){
        $output[4] = "NA";
    }
    
    # prognosis
    my @prognostics_value;
    for(my $i=56;$i<=72;$i++){
        if($.==1){
            push @prognostics_header, $l[$i];
        }else{
            push @prognostics_value, $l[$i];
        }
    }
    if($.>1){
        my $n=0;
        my $prognostics_value_out = "";
        foreach my $prognostics_value(@prognostics_value){
            if($prognostics_value =~ /(favourable .*|unfavourable .*)/){
                $prognostics_value = $1;
                if($prognostics_header[$n] =~ /Pathology prognostics - (.*)/){
                    if($prognostics_value_out eq ""){
                        $prognostics_value_out = "$1:".$prognostics_value;
                    }else{
                        $prognostics_value_out = $prognostics_value_out.",$1:".$prognostics_value;
                    }
                }            
            }
            $n++;
        }
        if($prognostics_value_out ne ""){
            $output[2] = $prognostics_value_out;
        }
    }
    
    # tissue expression
    my @expression_value;
    for(my $i=73;$i<=133;$i++){
        if($.==1){
            push @expression_header, $l[$i];
        }else{
            push @expression_value, $l[$i];
        }
    }
    if($.>1){
        my $m=0;
        my $expression_out = "";
        foreach my $expression_value(@expression_value){
            if($expression_header[$m] =~ /(heart)/ or $expression_header[$m] =~ /(kidney)/ or $expression_header[$m] =~ /(liver)/ or $expression_header[$m] =~ /(lung)/ or $expression_header[$m] =~ /(spleen)/ or $expression_header[$m] =~ /(midbrain)/ or $expression_header[$m] =~ /(cerebral cortex)/ or $expression_header[$m] =~ /(cerebellum)/){
                if($expression_out eq ""){
                    if($expression_value eq ""){
                        $expression_out = "$1:NA";
                    }else{
                        $expression_out = "$1:".$expression_value;
                    }
                }else{
                    if($expression_value eq ""){
                        $expression_out = $expression_out.",$1:NA";
                    }else{
                        $expression_out = $expression_out.",$1:".$expression_value;
                    }
                }
            }
            $m++;
        }
        $output[3] = $expression_out;
    }
    
    
    # output
    if($.>1){
        if(exists $ensembl_id{$l[2]}){
            $atlas{$ensembl_id{$l[2]}} = join("\t",@output);
        }
    }

}
close FILE5;


# hemogene
my $file6 = "./database/hemogene.txt";
my %hemogene;
my @hemo_header = ('AML','B-cells','CD34+','CD34+CD45RA-','Ery-I','Ery-II','Ery-III','Ery-IV','Gran-I','Gran-II','Gran-III','Gran-IV','Gran-V','Granulocytes','LSC_High','LSC_Low','LSC_Medium','Monocytes','nBM','PB_CD34+cells','Pre-B-I','Pre-B-II','T-cells','WBC');
open (FILE6, "< $file6") or die "Can't open $file6: $!";
while (my $l =<FILE6>){
    
    chomp $l;
    my @l = split(/\t/, $l);
    my @m = split(/,/, $l[1]);
    if(exists $ensembl_id{$l[0]}){
        my @hemo_out;
        for(my $i=0;$i<=23;$i++){
            push @hemo_out, $hemo_header[$i].":".$m[$i];
        }
        $hemogene{$ensembl_id{$l[0]}} = join(",",@hemo_out);;
    }

}
close FILE6;



### Create outputs

my $output = $output_path.$file_name.".SPAT_annotation.txt";
open (OUT, "> $output") or die "Can't open $output: $!";
print (OUT "INPUT\tGENE_NAME\tSPAT_SCORE\tPUBMED_IDs\tRED_FLAG_COUNT\tRED_FLAG_TERMS\tSECRETED_FLAG_COUNT\tSECRETED_FLAG_TERMS\tMAIN_LOCATION_NEXTPROT\tMAIN_LOCATION_HPA\tAPPROVED_ANTIBODY(see https://antibodyregistry.org/)\tSIGNIFICANT_PROGNOSTIC_VALUE_HPA\tTISSUE_RNA_HPA(NX)\tTISSUE_PROT_HPA(level)\tHEMOGENE_RNA(TPM)\n");

my $output2 = $output_path.$file_name.".log.out";
open (OUT2, "> $output2") or die "Can't open $output2: $!";
my $file7 = $input;
open (FILE7, "< $file7") or die "Can't open $file7: $!";
my %out; my $i=0; my $j=0; my $info=50;

while (my $l =<FILE7>){

    $l =~ s/\r\n?/\n/g; 
    chomp $l;
        
        my $ndb = $l;
        
        # test gene in db
        if(exists $transl{$l}){
            
            my @scores;
            
            # best score option
            if($best_score eq "true"){
            
                foreach my $ndb(@{$transl{$l}}){
                    my @l = split(/\t/, $spat{$ndb});
                    push @scores, $l[0];
                }
                
                # find best score in the scores array
                my $idx = firstidx { $_ eq max(@scores) } @scores;

                # output the corresponding entry
                $ndb = ${$transl{$l}}[$idx];
                
                if(exists $spat{$ndb}){
                        
                    $out{$ndb}=$l."\t".$gene_name{$ndb}."\t".$spat{$ndb};
                        
                    if(exists $atlas{$ndb}){
                        $out{$ndb}=$out{$ndb}."\t".$atlas{$ndb};
                    }else{
                        $out{$ndb}=$out{$ndb}."\tNA\tNA\tNA\tNA\tNA";
                    }
                        
                    if(exists $hemogene{$ndb}){
                        $out{$ndb}=$out{$ndb}."\t".$hemogene{$ndb};
                    }else{
                        $out{$ndb}=$out{$ndb}."\tNA";
                    }
                        
                    $i++;
                    
                }else{
                    $j++;
                }
                
            
            }else{
            
        
                foreach my $ndb(@{$transl{$l}}){
                    
                    if(exists $spat{$ndb}){
                        
                        $out{$ndb}=$l."\t".$gene_name{$ndb}."\t".$spat{$ndb};
                        
                        if(exists $atlas{$ndb}){
                            $out{$ndb}=$out{$ndb}."\t".$atlas{$ndb};
                        }else{
                            $out{$ndb}=$out{$ndb}."\tNA\tNA\tNA\tNA\tNA";
                        }
                        
                        if(exists $hemogene{$ndb}){
                            $out{$ndb}=$out{$ndb}."\t".$hemogene{$ndb};
                        }else{
                            $out{$ndb}=$out{$ndb}."\tNA";
                        }
                        
                        $i++;
                        
                    }else{
                        $j++;
                    }
                    
                }
                
            }
            
        }else{
            $j++;
        }
        
        # print count of annotated protein
        if($. == 1){print "\nSTARTING ANNOTATION:";}
        if($. == $info){
            print "\n$info proteins annotated...";
            $info = $info + 50;
        }

}
close FILE7;

foreach my $key(sort keys %out){

    my @l = split(/\t/, $out{$key});
    if($l[2] >= $score_threshold){
        print (OUT $out{$key}."\n");
    }
}
close OUT;



### Print message to log file

my $tot = $i+$j;
my $perc_anno = sprintf("%.2f",($i/($tot))*100);
my $perc_missing = sprintf("%.2f",($j/($tot))*100);
print (OUT2 "missing\t$j
total\t$tot 
annotated(%)\t$perc_anno\n");


### Create R file (and figure)

my $output3 = $output_path.$file_name.".SPAT_distribution.R";
open (OUT3, "> $output3") or die "Can't open $output3: $!";
print (OUT3 'suppressMessages(library(ggplot2))'."\n".'surf = read.table("'.$output.'", sep="\t", header = TRUE)'."\n".'g = ggplot(surf, aes(SPAT_SCORE)) + geom_density(adjust=1.2)'."\n".'y_max = ggplot_build(g)$layout$panel_scales_y[[1]]$range$range[2]'."\n".'png("'.$output_path.$file_name.'.SPAT_distribution.png", width=500, height=500, res=200)'."\n".'ggplot(surf, aes(SPAT_SCORE)) + geom_density(color="black", lwd=1, fill="grey", alpha=0.6, adjust=1.2) + geom_vline(aes(xintercept=mean(SPAT_SCORE)), color="#69b3a2", linetype="dashed", size=.7) + labs(x = "SPAT score") + scale_color_hue(name="Type") + theme_classic() + theme(plot.title = element_text(size=10)) + labs(title="SCORE DISTRIBUTION") + geom_label(label=paste("mean=",format(round(mean(surf$SPAT_SCORE), 2), nsmall = 2),sep=""), x = mean(surf$SPAT_SCORE), y = y_max, size = 1.8, label.padding = unit(0.2, "lines"), color = "black", fill="#69b3a2")'."\n".'dev.off()'."\n");  
close OUT3;

if($graph eq "true"){
    `Rscript $output_path/$file_name.SPAT_distribution.R`;
}


### Print final message

print "\n...\nSURFACE PROTEIN ANNOTATION COMPLETED\n\n";



exit 0;
       
