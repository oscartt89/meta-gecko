#!/usr/local/bin/perl
####Programa para obtener los locus_tag a partir de un genbank.

if (!$ARGV[0]){
        print "Introduzca el nombre del su fichero genbank con extension: ";
		$name=<STDIN>;
		chomp $name;
}
else{
	$name=$ARGV[0];
	}
	open (gbk,"$name") || die "Error: problem opening codigos\n";
	open (results, ">$name-parseado") || die "Error: problem creating results\n";

while (<gbk>){
	chomp $_;

	next unless (($_=~/^\s+\/locus_tag/) || ($_=~/^\s+\/db_xref="GI:/));
	if ($_=~/^\s+\/locus_tag/){
		($kk,$idtair)=split(/\=/,$_);
	#	print results "$idtair\t";
		}elsif($_=~/^\s+\/db_xref/){
			($qq,$idgi)=split(/\:/,$_);
			print results "$idtair\tgi|$idgi\n";
			}
}
close (gbk);
close (results);


				     


