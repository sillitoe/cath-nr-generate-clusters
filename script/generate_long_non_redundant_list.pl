#!/usr/bin/env perl

use strict;
use warnings;

use local::lib;
use Carp                       qw/ confess        /;
use Data::Dumper;
use English                    qw/ -no_match_vars /;
use feature                    qw/ say            /;
use List::MoreUtils            qw/ first_index    /;
use List::Util                 qw/ max min        /;
use Path::Class                qw/ dir file       /;
use Getopt::Long;

my $seq_id_cutoff  = 40.0;
my $overlap_cutoff = 60.0;
my $out_file;

my $program_name = file( $PROGRAM_NAME )->basename;

my $USAGE = <<"_USAGE";

usage: $program_name [options] <S100_domain_list> <blast_results_directory>

options:

  -o|--out   <file>          Output NR file [default: nr_list_s<SEQ>_o<OV>.txt]
  --seq      Num[20-100]     Specify sequence id cutoff [default: 40.0] 
  --overlap  Num[20-100]     Specify overlap cutoff [default: 60.0] 

_USAGE

GetOptions(
	'seq=s'     => \$seq_id_cutoff, 
	'overlap=s' => \$overlap_cutoff, 
	'o|out=s'   => \$out_file,
) or die "! Error parsing command line params\n$USAGE";

# Check argument
if (scalar(@ARGV) != 2) {
	die $USAGE;
}

if ( $seq_id_cutoff !~ /^[0-9\.]+$/ || $seq_id_cutoff < 20 || $seq_id_cutoff > 100 ) {
	die "! Error: expected --seq to be a number between 20 and 100 (not: $seq_id_cutoff)\n";
}

if ( $overlap_cutoff !~ /^[0-9\.]+$/ || $overlap_cutoff < 20 || $overlap_cutoff > 100 ) {
	die "! Error: expected --ov to be a number between 20 and 100 (not: $overlap_cutoff)\n";
}

$out_file = file( $out_file || sprintf( "nr_list_s%02d_o%02d.txt", $seq_id_cutoff, $overlap_cutoff) );

# Check arguments
my $s100_list_file = file( $ARGV[0] );
if ( ! -s $s100_list_file ) {
	confess "No such file \"$s100_list_file\"";
}
my $results_dir = dir( $ARGV[1] );
if ( ! -d $results_dir ) {
	confess "No such directory \"$results_dir\"";
}

# Parse S100s
my $sorted_s100_ids = parse_s100_ids_from_s100_domain_list($s100_list_file);

# Build the main data structure: a hash of hash refs,
#                                which stores (both symetric) links between pairs of domains as:
#                                    $data{$domain_id1}->{$domain_id2} = 1;
#                                    $data{$domain_id2}->{$domain_id1} = 1;
my %data = map { ( $ARG, {} ) } sort(@$sorted_s100_ids);
parse_blast_results($results_dir, \%data, $seq_id_cutoff, $overlap_cutoff);

# Choose the ids in the NR list 
my $nr_list = choose_nr_ids($sorted_s100_ids);

# Write the list of non-redundant doamins to a file (here "/tmp/nr_list.txt")
$out_file->spew(join("\n", sort(@$nr_list) )."\n");

# Grab the domain IDs from a CathDomainList.S100 file
sub parse_s100_ids_from_s100_domain_list {
	my $s100_list_file = shift;
	warn localtime(time())." : Reading S100 list from S100_domain_list_file $s100_list_file\n";
	my @s100_lines = $s100_list_file->slurp();
	while (chomp(@s100_lines)) {}
	my @s100_ids = map { (split( /\s+/, $ARG))[0] } sort(@s100_lines);
	@s100_ids = sort(@s100_ids);
	return \@s100_ids;
}

sub parse_blast_results {
	my $blast_results_dir = shift;
	my $data              = shift;
	my $seq_id_cutoff     = shift;
	my $overlap_cutoff    = shift;
	
	warn localtime(time())." : Parsing blast results from directory $blast_results_dir\n";
	
	# Loop over files in the directory
	my $results_file_ctr = 1;
	while (my $results_file = $results_dir->next) {
		next unless -f $results_file;
		
		# Grab the lines from the file
		my @results_lines = $results_file->slurp();
		while (chomp(@results_lines)) {}
		
		# Loop over the lines
		foreach my $results_line (@results_lines) {
			# Split the line on whitespace, check the number fo fields and assign to variables
			my @results_line_parts = split(/\s/, $results_line);
			if (scalar (@results_line_parts) != 6) {
				confess "Format of $results_file is invalid";
			}
			my ($id1, $id2, $seq_id, $aligned_length, $length1, $length2) = @results_line_parts;
			
			# Process the IDs and overlap
			$id1 = strip_domain_id($id1);
			$id2 = strip_domain_id($id2);
			my $overlap_over_longer = 100.0 * $aligned_length / max($length1, $length2);
			
			# Skip this entry if it doesn't meet the criteria:
			#  * not a self-hit
			#  * sequence identity  >= $seq_id_cutoff
			#  * percentage overlap >= $overlap_cutoff
			if ($id1 eq $id2 || $seq_id < $seq_id_cutoff || $overlap_over_longer < $overlap_cutoff) {
				next;
			}
			
			# Store the connection between these two results
			$data->{$id1}->{$id2} = 1;
			$data->{$id2}->{$id1} = 1;
		}
		
		# Output some progress every 1000 files
		if ($results_file_ctr % 1000 == 0) {
			warn localtime(time())." : Reading blast results file $results_file_ctr ($results_file)\n";
		}
		++$results_file_ctr;
# 		last if ($results_file_ctr > 1000);
	}
}

# Convert from ID like "cath|4_0_0|1xfcB02/18-233" to "1xfcB02"
sub strip_domain_id {
	my $id = shift;
	$id =~ s/^cath\|[\d_]+\|//g;
	$id =~ s/^cath\|current\|//g;
	$id =~ s/\/[\w\-\(\)]+$//g;
	return $id;
}

sub choose_nr_ids {
	my $sorted_ids = shift;
	
	my @nr_list;
	my %count_of_id = map { ( $ARG, scalar(keys( %{ $data{$ARG} } ) ) ) } sort(keys(%data));
	
	# While there is still data left
	while (scalar(keys(%data))) {
		
		# Find the first ID with the fewest number of neighbours
		my $fewest_neighbours;
		my $id_with_fewest_neighbours;
		foreach my $id (@$sorted_ids) {
			my $num_neighbours = $count_of_id{$id};
			if (!defined($fewest_neighbours) || $num_neighbours < $fewest_neighbours) {
				$fewest_neighbours         = $num_neighbours;
				$id_with_fewest_neighbours = $id;
			}
			last if ($fewest_neighbours == 0);
		}
		my $id = $id_with_fewest_neighbours;
		
		# Grab the neighbours of this ID
		my @neighbours = sort(keys(%{$data{$id}}));
		
		# Output some information about this choice
		warn localtime(time())." : Next ID from list with ".scalar(keys(%data))." entries is $id with ".scalar(@neighbours)." neighbours\n";
		
		# Remove any connections to the neighbours of this ID
		foreach my $neighbour (@neighbours) {
			my @neighbours_of_neighbour = sort(keys(%{$data{$neighbour}}));
			foreach my $neighbour_of_neighbour (@neighbours_of_neighbour) {
				delete $data{$neighbour_of_neighbour}->{$neighbour};
				--($count_of_id{$neighbour_of_neighbour});
			}
		}
		# Remove this ID and its neighbours
		foreach my $id_or_neighbour (@neighbours, $id) {
			delete $data{$id_or_neighbour};
			delete $count_of_id{$id_or_neighbour};
			my $index_in_sorted_ids = first_index {$ARG eq $id_or_neighbour} @$sorted_ids;
			if ( $index_in_sorted_ids == -1 ) {
				confess "Could not find neighbour/id \"$id_or_neighbour\" in list of sorted ids (size " . scalar( @$sorted_ids ) . ") whilst processing id \"$id\"";
			}
			splice(@$sorted_ids, $index_in_sorted_ids, 1);
		}
		
		# Add the ID to the non-redundant list
		push @nr_list, $id;
	}
	
	@nr_list = sort(@nr_list);
	
	return \@nr_list;
}
