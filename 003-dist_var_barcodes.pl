$fname = $ARGV[0];
$true_barcode_file = $ARGV[1];


open FILE, $fname or die "can't open $fname\n";
chomp(@bc = <FILE>);
close FILE;


open FILE, $true_barcode_file or die;
chomp(@true_bc = <FILE>);
close FILE;


sub hd{ length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

$out = $fname . ".out";
open OUT, ">$out" or die;
print OUT "idx_true_bc_1\tdist1\tidx_true_bc_2\tdist2\n";

foreach $bc (@bc){
	
	@h = map{hd($_, $bc)} @true_bc;
	
	# we get the index for the closest
	$idxMin = 0;
	$h[$idxMin] < $h[$_] or $idxMin = $_ for 1 .. $#h;
	$ham1 = $h[$idxMin];
	
	# now we want the second closest
	$h[$idxMin] = 1000; # this ensures that the minimum value is not minimum again
	$idxMin2 = 0;
	$h[$idxMin2] < $h[$_] or $idxMin2 = $_ for 1 .. $#h;
	$ham2 = $h[$idxMin2];
	
	#convert indexes to R numbering
	$idxMin++;
	$idxMin2++;
	
	print OUT "$idxMin\t$ham1\t$idxMin2\t$ham2\n";
	
}
close OUT;

