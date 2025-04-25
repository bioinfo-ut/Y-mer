# #   #
#  # #     ##### ###  ##
#   #  ### # # # #-  #
#   #      # # # ### #

# REQUIED SOFTWARE
# PERL
# GENOMETESTER4 PACKAGE https://github.com/bioinfo-ut/GenomeTester4
# R

#REQUIRED R, PERL SCRIPTS AND SINGLE COPY K-MERS
# add_compare_depth_k_mers.pl    #FINDS INTERSEQTION BETWEEN MODEL AND DEPTH K-MERS AND EXCLUDES THEM
# distribution.pl                #CALCULATING TRESHOLD FOR SEQUENCING ERRORS AND MEDIAN COVERAGE
# MWS.R                          #MANN-WHITNEY TEST
# MODEL.R                        #MODEL CREATOR
# Y_NIPT_nimekiri.txt            #SINGLE COPY K-MERS IN CHRY (GRCH38)

# WHIS PIPELINE BUILD Y-MER MODEL FOR PREDICTING CHROMOSOME Y HAPLOGROUP (HG) INCLUDED IN MALES LISTING (males.txt)
# 5 STEPS ARE: CREATING FEMALES K-MERS UNION LIST, CREATING HG-S INTERSECTION LISTS, CREATING CHRY K-MERS LIST,
# FINDING INFORMATIVE K-MERS AND BUILDING THE MODEL


#
$gtester = "GenomeTester4/src";    #we use GenomeTester4 for k-mer manipulations
$working = "lists";                #SSD disk with best ReadWrite parameters
$dsR = "data";                     #R scripts, calculated k-mer counts
$lists = "lists";                  #k-mer binary list locations
$mens = "men.txt";                #one group of samples
$womens = "women.txt";            #second group of samples


# 1. FEMALE UNION 
# STEPS 
# A. MAKING FASTQ FILE, 
# B. MAKING LIST FILE,
# C. REMOVING SEQUENCING ERRORS BY EXCLUDING K-MERS WITH LOWER FREQUENCIES 
# D. ADDING FEMALES K-MERS TO UNION LIST

open SISSE, "$womens" or die;
while(<SISSE>){
   chomp;
   @tmp = split(/\t/);
   for($i = 1; $i < scalar(@tmp); $i++){
      system("cp ".$tmp[$i].".aligned.bam* $working"); # BAM FILE COPING, (WGET FROM WWW)
      system("samtools bam2fq $working/".$tmp[$i].".aligned.bam > $working/".$tmp[$i].".fastq");   #A
      system("rm $working/".$tmp[$i].".aligned.bam*");
      system("$gtester/glistmaker $working/".$tmp[$i].".fastq -w 25 -o $working/".$tmp[$i]."");    #B
      ### C EXCLUDING K-MERS WITH LOWER FREQUENCIES
      system("$gtester/glistquery $working/".$tmp[$i]."_25.list --distribution 100 |perl distribution.pl |head -1 > $lists/".$tmp[$i]."_25.txt");
      open JAOTUS, "$lists/".$tmp[$i]."_25.txt" or die;
      while(<JAOTUS>){
         chomp;
         @jaotus = split(/\t/);
         $katvus{$tmp[$i]} = $jaotus[1];
         system("$gtester/glistcompare $working/".$tmp[$i]."_25.list $working/".$tmp[$i]."_25.list -i -c $jaotus[0] -o $working/".$tmp[$i].""); #C
         system("rm $working/".$tmp[$i]."_25.list");
      }
      close JAOTUS;
      ####

      $k = $i -1;
      system("$gtester/glistcompare $working/".$tmp[$i]."_25_intrsec.list $working/".$tmp[$i]."_25_intrsec.list -u -o $working/".$tmp[0]."_".$i."") if $i == 1; #D
      system("$gtester/glistcompare $working/".$tmp[0]."_".$k."_25_union.list $working/".$tmp[$i]."_25_intrsec.list -u -o $working/".$tmp[0]."_".$i."") if $i > 1;  #D
      system("rm $working/".$tmp[0]."_".$k."_25_union.list");
      system("mv $working/".$tmp[$i]."_25_intrsec.list /gpfs/space/GI/GV/Projects/Y-mer/pb/females");
   }
   $i = $i -1;

   ### FEMALES UNION LIST
   system("mv $working/".$tmp[0]."_".$i."_25_union.list $working/".$tmp[0]."_25_union.list");
   system("cp $working/".$tmp[0]."_25_union.list /gpfs/space/GI/GV/Projects/Y-mer/pb/females/");

}
close SISSE;

# 2. MALE PER HGs, INTERSECTION
# STEPS 
# A. MAKING FASTQ FILE, 
# B. MAKING LIST FILE 
# C. REMOVING SEQUENCING ERRORS BY EXCLUDING K-MERS WITH LOWER FREQUENCIES 
# D. ADDING MALES K-MERS TO HG INTERSECTION LIST

open SISSE, "$mens" or die;
while(<SISSE>){
   chomp;
   @tmp = split(/\t/);
   push @grupid, $tmp[0];
   for($i = 1; $i < scalar(@tmp); $i++){
      system("cp ".$tmp[$i].".aligned.bam* $working");
      system("samtools bam2fq $working/".$tmp[$i].".aligned.bam > $working/".$tmp[$i].".fastq"); #A
      system("rm $working/".$tmp[$i].".aligned.bam*");
      system("$gtester/glistmaker $working/".$tmp[$i].".fastq -w 25 -o $working/".$tmp[$i]."");  #B
      ### C EXCLUDING K-MERS WITH LOWER FREQUENCIES
      system("$gtester/glistquery $working/".$tmp[$i]."_25.list --distribution 100 |perl distribution.pl |head -1 > $lists/".$tmp[$i]."_25.txt");
      open JAOTUS, "$lists/".$tmp[$i]."_25.txt" or die;
      while(<JAOTUS>){
         chomp;
         @jaotus = split(/\t/);
         $katvus{$tmp[$i]} = $jaotus[1];
         system("$gtester/glistcompare $working/".$tmp[$i]."_25.list $working/".$tmp[$i]."_25.list -i -c $jaotus[0] -o $working/".$tmp[$i].""); #C
         system("rm $working/".$tmp[$i]."_25.list");
      }
      close JAOTUS;
      ####
      $k = $i -1;
      system("$gtester/glistcompare $working/".$tmp[$i]."_25_intrsec.list $working/".$tmp[$i]."_25_intrsec.list -i -o $working/".$tmp[0]."_".$i."") if $i == 1; #D
      system("$gtester/glistcompare $working/".$tmp[0]."_".$k."_25_intrsec.list $working/".$tmp[$i]."_25_intrsec.list -i -o $working/".$tmp[0]."_".$i."") if $i > 1; #D
      system("rm $working/".$tmp[0]."_".$k."_25_intrsec.list");
      system("cp $working/".$tmp[$i]."_25_intrsec.list /gpfs/space/GI/GV/Projects/Y-mer/pb/males");
   }
   $i = $i -1;

   ### MALES HGs INTERSECTION LIST
   system("mv $working/".$tmp[0]."_".$i."_25_intrsec.list $working/".$tmp[0]."_25_intrsec.list");
   ### MALES HGs INTERSECTION LIST WITHOUT FEMALES K-MERS IE. CHRY SPECIFIC HGs K-MERs LISTS
   system("$gtester/glistcompare $working/".$tmp[0]."_25_intrsec.list $working/fem_25_union.list -dd -o $working/".$tmp[0]."");
   system("rm $working/".$tmp[0]."_25_0_diff2.list");
}
close SISSE;



# 3. MALES HGs UNION
# ALL CHRY SPECIFIC HGs K-MERs LISTS ARE JOINED TO MALES MODEL HSs UNIOIN
$i = 0;
open SISSE, "$mens" or die;
while(<SISSE>){
   chomp;
   @tmp = split(/\t/);
   $i++;
   $k = $i -1;
   system("$gtester/glistcompare $working/".$tmp[0]."_25_0_diff1.list $working/".$tmp[0]."_25_0_diff1.list -u -o $working/MENS_".$i."") if $i == 1;
   system("$gtester/glistcompare $working/MENS_".$k."_25_union.list $working/".$tmp[0]."_25_0_diff1.list -u -o $working/MENS_".$i."") if $i > 1;
   system("rm $working/MENS_".$k."_25_union.list");
}

system("mv $working/MENS_".$i."_25_union.list $working/HGs_FEMOUT_25_union.list");
close SISSE;


## MALES MODEL HSs UNIOIN IN TXT FILE WITHOUT FREQUENCIES FOR USING GLISTQUERY
system("$gtester/glistquery $working/HGs_FEMOUT_25_union.list | cut -f 1 > $lists/HGs_FEMOUT_25_union.txt");




# 4.1. COUNTS PER MALE FOR FMCENTRO
# USING GLISTQUERY ALL K-MERS FREQUENCIES ARE QUERIED FOR ALL MALES SAMPLES AND NORMALIZED WITH SEQUENCING DEPTH CALCULATED WITH HELP OF JAOTUS.PL

open SISSE, "$mens" or die;
while(<SISSE>){
   chomp;
   @tmp = split(/\t/);
   for($i = 1; $i < scalar(@tmp); $i++){
     system("$gtester/glistquery $working/".$tmp[$i]."_25_intrsec.list -f $lists/HGs_FEMOUT_25_union.txt |cut -f 2 > $dsR/".$tmp[$i].".counts1");
     open MAN,"$dsR/".$tmp[$i].".counts1" or die;
     open MAN2,">$dsR/".$tmp[$i].".counts1_norm" or die;
     while(<MAN>){
        chomp;
        $freq = int($_/$katvus{$tmp[$i]}*10)/10; #NORMALIZING
        print MAN2 "$freq\n";
     }
     close MAN;
     close MAN2;
   }
}
close SISSE;


# 4.2. HELPING FILES FOR MWS INPUT

open K_MER,"$lists/HGs_FEMOUT_25_union.txt" or die;
open K_MER_2, ">$lists/HGs_FEMOUT_25_union_tulp1.txt" or die;
open K_MER_TAB, ">$lists/HGs_FEMOUT_25_union_tab.txt" or die;
print K_MER_2 "\n\n";
print K_MER_TAB "\n\n";
while(<K_MER>){
  chomp;
  print K_MER_TAB "\n";
  print K_MER_2 "$_\n";
}
close K_MER_TAB;
close K_MER_2;
close K_MER;

system("tail -n +2 $lists/HGs_FEMOUT_25_union_tab.txt > $lists/HGs_FEMOUT_25_union_tab_2.txt");

# 4.3. FILES FOR MWS AND MWS
# CALCULATES MANN-WHITNEY TEST FOR EVERY K-MER BASED ON AVERAGE FREQUENCY IN HG AND OUT OF HG SAMPLES

open NIMED, "$mens" or die;
open PAIS_1, ">FM_inid.txt" or die; 
while(<NIMED>){
   chomp;
   @tmp = split(/\t/);
   for($i = 1; $i < scalar(@tmp); $i++){
      print PAIS_1 "\t$tmp[$i]";
      $kood{$tmp[$i]} = $tmp[0];
      $paste .= " $dsR/";
      $paste .= "$tmp[$i]";
      $paste .= ".counts1_norm";
   }
}
print PAIS_1 "\n";

system("paste $lists/HGs_FEMOUT_25_union_tab_2.txt ".$paste." |head -n -1 > $lists/HGs_FEMOUT_25_union_tabel.txt");
close NIMED;

## CYCLE FOR EVERY HG
foreach $grupp (@grupid){
    open PAIS_2, ">".$grupp."_inid.txt" or die;
    open NIMED, "$mens" or die;
    while(<NIMED>){
       chomp;
       @tmp = split(/\t/);
       for($i = 1; $i < scalar(@tmp); $i++){
          print PAIS_2 "\tM" if $kood{$tmp[$i]} eq $grupp;
          print PAIS_2 "\tF" if $kood{$tmp[$i]} ne $grupp;
       }
    }
    print PAIS_2 "\n";
    system("cat FM_inid.txt ".$grupp."_inid.txt $lists/HGs_FEMOUT_25_union_tabel.txt > sisse.txt");
    system("paste $lists/HGs_FEMOUT_25_union_tulp1.txt $lists/HGs_FEMOUT_25_union_tab.txt $lists/HGs_FEMOUT_25_union_tab.txt $lists/HGs_FEMOUT_25_union_tab.txt $lists/HGs_FEMOUT_25_union_tab.txt $lists/HGs_FEMOUT_25_union_tab.txt sisse.txt > ".$grupp."_sisse.txt");
    system("Rscript MWS.R ".$grupp."_sisse.txt ".$grupp."_valja.txt");
    close NIMED;
    close PAIS_2;
}



# 4.4. MWS OUTPUT SORTING BY SPECIFITY IN HG


foreach $grupp (@grupid){
    system("cut -d ' ' -f 1-4 ".$grupp."_valja.txt | sort -k 2 -n -r > ".$grupp."_sort.txt");
}

# 4.5. SELECTING OF K-MER IN STEP OF 10000


@arvud = ('50000');
foreach $arv (@arvud){
   system("rm k-merid_".$arv.".txt");
   foreach $grupp (@grupid){
      system("cut -c 1-25 ".$grupp."_sort.txt |head -".$arv." >> k-merid_".$arv.".txt");
      system("cut -c 1-25 ".$grupp."_sort.txt |tail -".$arv." >> k-merid_".$arv.".txt");
   }
   
system("sort -u k-merid_".$arv.".txt > k-merid_unic_".$arv.".txt");
## COMPARISION OF MODEL AND SEQUENCING DEPTH K-MERS IN SETS AND EXCLUDING PRESENTED IN BOTH SETS, ADDING "M" AND "N" MARKING FOR DIFFERENTIATE MODEL AND DEPTH K-MERS
system("perl add_compare_depth_k_mers.pl k-merid_unic_".$arv.".txt > k-merid_unic_".$arv."_NIPT.db");
## MODEL K-MERS IN TXT FILE WITHOUT FREQUENCIES FOR USING GLISTQUERY
system("cut -f 3 k-merid_unic_".$arv."_NIPT.db > k-merid_unic_".$arv."_NIPT.txt");

}

# CREATING MODEL
# 5.1. COUNTS FOR MODEL INPUT K-MERS

open SISSE, "$mens" or die;
while(<SISSE>){
   chomp;
   @tmp = split(/\t/);
   foreach $arv (@arvud){
      for($i = 1; $i < scalar(@tmp); $i++){
         system("$gtester/glistquery $working/".$tmp[$i]."_25_intrsec.list -f k-merid_unic_".$arv."_NIPT.txt |cut -f 2 > $dsR/".$tmp[$i]."_".$arv.".counts2");
      }
   }
}
close SISSE;

## 5.2. CREATING MODEL INPUT FILE AND CALCULATING MODEL

## CYCLE FOR DIFFERENT SET OF K-MERS (STEP OF 10000) IF NEEDED
foreach $arv (@arvud){
   open KMER, "k-merid_unic_".$arv."_NIPT.db" or die;
   $rida = 2;
   while(<KMER>){
      chomp;
      @tmp1 = split(/\t/);
      $tabel[$rida][0] = $tmp1[0];
      $rida++;
   }
   close KMER;

 
   open NIMED, "$mens" or die;
   open VALJA, ">input_nimekiri_".$arv.".txt" or die;
   $tulp = 1;
   while(<NIMED>){
      chomp;
      @tmp = split(/\t/);
      for($i = 1; $i < scalar(@tmp); $i++){
         $tabel[0][$tulp] = $tmp[$i];
         $tabel[1][$tulp] = $tmp[0];
         $depth = $katvus{$tmp[$i]};
         $rida = 2;
   
         open IN, "$dsR/".$tmp[$i]."_".$arv.".counts2" or die;
         while(<IN>){
            chomp;
            $tabel[$rida][$tulp] = $_;
            $rida++;
         }
         close IN;
         $tulp++;
      }
   }
   close NIMED;
   #mudelisse algab $t = 1, callimiseks $t = 2

   #print VALJA "$tabel[0][1]";
   for($t = 1; $t < $tulp ; $t++){
      print VALJA "\t$tabel[0][$t]";
   }
   print VALJA "\n";
   for($t = 1; $t < $tulp ; $t++){
      print VALJA "\t$tabel[1][$t]";    
   }
   print VALJA "\n"; 
         
   for($r = 2; $r < $rida; $r++){
      print VALJA "$tabel[$r][0]"; #siin ka 0 vs.1
      for($t = 1; $t < $tulp ; $t++){
         print VALJA "\t$tabel[$r][$t]";
      }
      print VALJA "\n";
   }
   close VALJA;

### CREATING MODEL
   system("Rscript MODEL.R input_nimekiri_".$arv.".txt mudel_nimekiri_".$arv.".Rdata");
}
