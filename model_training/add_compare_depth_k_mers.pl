# Informative k-mers for model "m" and for depth estimation "n" are marked if exclude if k-mer is in both set

#informative k-mers from all HG-s (may be multiplicated)
open AUGUST, "$ARGV[0]" or die;
while(<AUGUST>){
   chomp;
   $august{$_} = "m";
}
close AUGUST;

# single copy k-mers presented only in Y chromosome
open NIPT, "Y_NIPT_nimekiri.txt" or die;
while(<NIPT>){
   chomp;
   $nipt{$_} = "n";
}
close NIPT;


foreach $r (keys %august){
   print"m$r\t1\t$r\n" unless exists $nipt{$r} || $nipt{&rev_compl($r)};
}
foreach $p (keys %nipt){
   print"n$p\t1\t$p\n" unless exists $august{$p} || $august{&rev_compl($p)};
}

sub rev_compl {
   my $seqref = shift;
   my $revcomp = reverse $seqref;
   $revcomp =~ tr [AaCcGgTt] [TtGgCcAa];
   return $revcomp;
}
                