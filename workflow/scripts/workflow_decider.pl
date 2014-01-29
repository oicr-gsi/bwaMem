use strict;
use XML::DOM;
use Data::Dumper;

 my $parser = new XML::DOM::Parser;
 #my $doc = $parser->parsefile ("https://cghub.ucsc.edu/cghub/metadata/analysisDetail?participant_id=3f70c3e3-0131-466f-92aa-0a63ab3d4258");
 # FIXME: why is this file truncated?
 #system("lwp-download 'https://cghub.ucsc.edu/cghub/metadata/analysisDetail?study=TCGA_MUT_BENCHMARK_4&state=live' data.xml");
 #my $doc = $parser->parsefile ('https://cghub.ucsc.edu/cghub/metadata/analysisDetail?study=TCGA_MUT_BENCHMARK_4&state=live');
 #######system("cgquery -o data.xml \"study=TCGA_MUT_BENCHMARK_4&state=live\"");
 my $doc = $parser->parsefile("data.xml");

 # print all HREF attributes of all CODEBASE elements
 my $nodes = $doc->getElementsByTagName ("Result");
 my $n = $nodes->getLength;

 for (my $i = 0; $i < $n; $i++)
 {
     my $node = $nodes->item ($i);
     #$node->getElementsByTagName('analysis_full_uri')->item(0)->getAttributeNode('errors')->getFirstChild->getNodeValue;
     print $node->getElementsByTagName('analysis_full_uri')->item(0)->getFirstChild->getNodeValue;
     my $aurl = $node->getElementsByTagName('analysis_full_uri')->item(0)->getFirstChild->getNodeValue;
     ######system("wget -q -O data2.xml $aurl");
     my $adoc = $parser->parsefile ("data2.xml");

     #my $hrefs = $node->getElementByTagName ("library_strategy");
     #my $libCnt = $hrefs->getLenght;
     print "\n";

     my $libName = $adoc->getElementsByTagName('LIBRARY_NAME')->item(0)->getFirstChild->getNodeValue;
     my $libStrategy = $adoc->getElementsByTagName('LIBRARY_STRATEGY')->item(0)->getFirstChild->getNodeValue;
     my $libSource = $adoc->getElementsByTagName('LIBRARY_SOURCE')->item(0)->getFirstChild->getNodeValue;
     print "$libName $libStrategy $libSource\n";

     exit;
 }

 # Print doc file
 #$doc->printToFile ("out.xml");

 # Print to string
 #print $doc->toString;

 # Avoid memory leaks - cleanup circular references for garbage collection
 $doc->dispose;


