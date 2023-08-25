use strict;
use warnings;
use Cwd;

my $base = getcwd;
my @buckets = ("output-run2-protein");#,"output-run2-protein");#,"output-cgenff-water","output-cgenff-protein","output-openff-v2-water","output-openff-v2-protein");
#@buckets = ("output-gaff2-sh-water");#,"output-cgenff-protein","output-openff-v2-protein");
#my $bucket = $ARGV[0];

# eg s3://output-gaff2-sh-water/syk/edge_CHEMBL3265023_CHEMBL3265022/stateB/run3/metadata_1624888020

foreach my $bucket(@buckets)
{
	system("aws s3 ls s3://$bucket/ --recursive  > raw_bucket_$bucket");
	print("$bucket: Rsync metadata\n");
	system("mkdir $base/$bucket");
	system("grep metadata raw_bucket_$bucket > raw_metadata_$bucket");
        do_sync( "raw_metadata_$bucket", $bucket );

	system("rm raw_bucket_$bucket raw_metadata_$bucket");
}

sub do_sync
{
	open(FOO,$_[0]);
	my @cont = <FOO>;
	close(FOO);
	my $bucket = $_[1];

	foreach my $line(@cont)
	{
#		print("$line\n");
                chomp($line);
		my @arr = split(/[\s\t]{1,}/,$line);
		my @folders = split(/\//,$arr[-1]);

		my $path = "$bucket";		
		for(my $i=0; $i<=3; $i++)
		{
			$path .= "/$folders[$i]";
			system("mkdir $path");
		}
#		system("aws s3 sync s3://$path/ ./$path/. --exclude \"*\" --include \"metadata*\"\n");
#		system("aws s3 sync s3://$path/ ./$path/. --exclude \"*\" --include \"env*\"\n");
#		system("aws s3 sync s3://$path/ ./$path/. --exclude \"*\" --include \"dhdl79.xvg\"\n");
		system("aws s3 sync s3://$path/ ./$path/. --exclude \"*\" --include \"work.dat\"\n");
#		system("aws s3 sync s3://$path/ ./$path/. --exclude \"*\" --include \"md*log\"\n");
#		system("aws s3 sync s3://$path// ./$path/. --exclude \"*\" --include \"md*log\"\n");
			
	}
}


