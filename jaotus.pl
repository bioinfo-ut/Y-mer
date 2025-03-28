$rida = 0;
while(<>){
   chomp;
   @tmp = split(/\t/);
   $jrk[$tmp[0]] = $tmp[1];
   $jrk++;
}
for($i = 4; $i < $jrk; $i++){
   print "$i" if ($jrk[$i] <  $jrk[$i-1] && $jrk[$i] <  $jrk[$i+1]);
   print "\t$i\n" if ($jrk[$i] >  $jrk[$i-1] && $jrk[$i] >  $jrk[$i+1]);
}
   