open AUGUST, "$ARGV[0]" or die;
while(<AUGUST>){
   chomp;
   $august{$_} = "m";
}
close AUGUST;

open NIPT, "/storage7/inimene/CHM13/t2t-chm13-v2.0/INDEX_25/chrY/uuesti_august22/NIPT_v7/Y_NIPT_nimekiri.txt" or die;
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
                