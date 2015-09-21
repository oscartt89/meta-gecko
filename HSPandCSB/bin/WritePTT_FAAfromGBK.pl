#!/usr/bin/perl -w

my $gbk_file = $ARGV[0];
my $org_name = $ARGV[1];

if ($#ARGV < 1) {
  print "Usage: WritePTT_FAAfromGBK.pl  <gbk_file_name> <genome_name>\n";
  print "     (example: WritePTT_FAAfromGBK.pl NC_007332.gbk Mycoplasma_hyopneumoniae_7448)\n";
  exit(1);
}

$found_gene = "false";
$found_seq = "false";

#print "\nWriting files: $org_name.ptt and $org_name.faa\n";

open (PTT,">$org_name.ptt");

$org_name_spaces = $org_name;
$org_name_spaces =~ s/\_/ /g;

#open (FAA,">$org_name.faa");
my $gene = "-";
open(ARQ,"$gbk_file");
while (<ARQ>)
  {
    if (/ORGANISM\s+(.+)/)
      {
	$org = $1;
      }
    elsif (/^LOCUS\s+\S+\s+(\d+)\s+/) #LOCUS       CIAT899_chromosome   3837061 bp    DNA     linear       05-MAR-2012
      {
   	$genome_size = $1;
      }
    elsif (/CDS\s*(\d+)\.\.(\d+)/)
      {       
	$strand =  "f";
	$start = $1;
	$stop = $2;
	$gi = "";
	$found_gene = "true";	       
      }
    elsif (/CDS\s*complement\((\d+)\.\.(\d+)/)
      {       
	$strand =  "r";
	$start = $1;
	$stop = $2;
	$gi = "";
	$found_gene = "true";
      }
    elsif (/CDS\s+.*join\((\d+)\.\..+\.\.(\d+)/)
      {
        $strand =  "f";
        $start = $1;
        $stop = $2;
        $gi = "";
        $found_gene = "true";
      }
    elsif (($found_gene eq "true") && (/\/gene\=\"(.+)\"/))
      {       
	$gene = $1;
	if ($gene =~ /hypothetical/)
	  {
	    $gene = "-";
	  }
      }      
    elsif (($found_gene eq "true") && (/\s+\/product=\"(.+)\"/))
      {
	$product = $1;
      }
    elsif (($found_gene eq "true") && (/\s+\/protein_id=\"(.+)\"/))
      {
	#$syn = $1;
	$found_gene = "false";
	#$size = ($stop - $start - 2)/3; 
	#print NEW "\t$start..$stop\t$strand\t$size\t$syn\t$gene\t$syn\t\-\t\-\t$product\n";
	#$gene = "-";
      }
    elsif (/translation=\"(\w+)\"/)    #/translation="MDFTEEFSQIYKKEKCKIDTKSQFWMYMKQFNKKLKYN"
      {
	$seq = $1;
	$found_seq = "false";
	$gi = $syn if ((not defined $gi) || ($gi eq ""));
#	print FAA ">gi|$gi|ref|$syn| $product [$org_name_spaces]\n";
#	print FAA "$seq\n";
	$size = length($seq);
	#print PTT "$start..$stop\t$strand\t$size\t$gi\t$gene\t$syn\t\-\t\-\t$product\n";
	$ptt_line .= "$start..$stop\t$strand\t$size\t$gi\t$gene\t$syn\t\-\t\-\t$product\t\n";
	$cont_genes++;
	$gene = "-";
      }
    elsif (/translation=\"(\w+)/)   
      {
	$seq = $1;
	$found_seq = "true";      
      }
    #elsif (($found_seq eq "true") && (/\s+([A-Z]+)/))    
    elsif ((/\s+([A-Z]*)\"/) && ($found_seq eq "true"))
      {	
	$seq .= $1;
	$found_seq = "false";
	$gi = $syn if ((not defined $gi) || ($gi eq ""));
	#antes print FAA ">gi|$gi|ref|$syn|$gene|$product [$org]\n";
	#antes print FAA ">gi|$gi|ref|$syn| $product [$org_name_spaces]\n";
#	print FAA ">gi|$gi|ref|$syn| $product [$org]\n";
#	print FAA "$seq\n";
        $size = length($seq);	
	#print PTT "$start..$stop\t$strand\t$size\t$gi\t$gene\t$syn\t\-\t\-\t$product\n";
	$ptt_line .= "$start\t$stop\t$strand\t$size\t$gi\t$gene\t$syn\t\-\t\-\t$product\t\n";
	$cont_genes++;
	$gene = "-";
      }
    elsif ( (/\s+([A-Z]+)/) && ($found_seq eq "true"))
      {
	$seq .= $1;	
      }
    elsif (/\s+\/db_xref=\"GI\:(.+)\"/)
      {
	$gi = $1;
	$syn = $gi if ((not defined $syn) || ($syn eq ""));	
	#$found_seq = "false";
	#$gene = "-";
      }
    elsif (/\s+\/locus_tag=\"(.+)\"/)
      {
	$syn = $1;
	#$found_seq = "false";
	#$gene = "-";
      }	 
  }
print PTT "$org, complete genome - 1..$genome_size\n";
#print PTT "$org_name_spaces, complete genome - 1..$genome_size\n";
print PTT "$cont_genes proteins\n";
print PTT "Location         Strand  Length  PID     Gene    Synonym Code    COG     Product\n";
print PTT $ptt_line;

#close FAA;
close PTT;
