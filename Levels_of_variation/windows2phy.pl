#!/usr/bin/perl -w

use strict;
use warnings;
use 5.010;
use Data::Dumper;

open ( my $in , "<" , $ARGV[ 0 ] ) or die "$!";

my $header = <$in>;
chomp $header;

my @headers = split ( /\t/ , $header );

shift @headers;
shift @headers;

my @pops;

for ( my $i = 0; $i < $#headers; $i += 2 ) {

	my ( $n , $fst ) = @headers[ $i , $i + 1 ];

	if ( $fst =~ m/FST\_(\S+)\_vs\_(\S+)/ ) {
	
		push @pops , ( $1 , $2 );

	}
	
}

my $min_fst_table = {};

my @chroms = ( "all" );

my $chrom_table = {};

while (<$in>) {

	chomp;
	
	if ( $_ =~ m/^(\S+)\t(\S+)$/ ) {
	
		next;
	
	}
	
	my @data = split ( /\t/ , $_ );
	
	my $chrom = shift @data;
	my $pos = shift @data;
	
	next unless @data;
	
	my $fst_table = {};
	
	my $keep = 0;
	
	for ( my $i = 0; $i < $#data; $i += 2 ) {
	
		my ( $n , $fst ) = @data[ $i , $i + 1 ];
	
		if ( $n =~ m/\d/ ) {
		
			$keep = 1;
			
			last;
		
		}
	
	}
	
	next unless $keep;

	say STDERR $_;
	
	for ( my $i = 0; $i < $#data; $i += 2 ) {
	
		my ( $n , $fst ) = @data[ $i , $i + 1 ];
	
		unless ( $n =~ m/\d/ ) {
		
			$n = 0;
			$fst = 0;
		
		}
	
		my ( $pop , $other_pop ) = @pops[ $i , $i + 1 ];
		
		$fst_table->{ $pop }->{ $other_pop } = $fst;
		$fst_table->{ $other_pop }->{ $pop } = $fst;
		
	}
	
	my @sorted_pops = keys %{ $fst_table };
	
	say scalar @sorted_pops;
	
	for ( my $i = 0; $i < @sorted_pops; $i++ ) {
	
		my $pop = $sorted_pops[ $i ];
	
		print "$pop       ";
	
		my $min_fst = 1.1;
		my $min_fst_pop = "";
	
		for ( my $j = 0; $j < @sorted_pops; $j++ ) {
	
			my $other_pop = $sorted_pops[ $j ];
	
			if ( $pop eq $other_pop ) {
			
				print " " , sprintf ( "%.4f" , 0 );
			
			}
			
			elsif ( not defined $fst_table->{ $pop }->{ $other_pop } ) {
			
				print " " , sprintf ( "%.4f" , 0 );
			
			}
			
			else {
	
				print " " , sprintf ( "%.4f" , $fst_table->{ $pop }->{ $other_pop } );
	
				if ( $fst_table->{ $pop }->{ $other_pop }  < $min_fst ) {
		
					$min_fst = $fst_table->{ $pop }->{ $other_pop } ;
					$min_fst_pop = $other_pop;
		
				}
	
			}
	
		}
		
		if ( $min_fst <= 1 ) {
		
			$min_fst_table->{ $pop }->{ 'all' }->{ 'group' }->{ $min_fst_pop }++;
			$min_fst_table->{ $pop }->{ 'all' }->{ 'n' }++;
			
			push @chroms , $chrom unless defined $chrom_table->{ $chrom };
			$chrom_table->{ $chrom } = 1;
			
			$min_fst_table->{ $pop }->{ $chrom }->{ 'group' }->{ $min_fst_pop }++;
			$min_fst_table->{ $pop }->{ $chrom }->{ 'n' }++;
			
		}

		say "";
	
	}
	
	say "";
	
	# exit;

}

foreach my $pop ( sort keys %{ $min_fst_table } ) {

	open ( my $out , ">" , $ARGV[ 0 ] . ".closest_group.${pop}.tsv" ) or die "$!";
	
	say $out "GROUP\tCHROM\tWINDOWS\t" , join ( "\t" , sort keys %{ $min_fst_table } );
	
	foreach my $chrom ( @chroms ) {
	
		my $n = $min_fst_table->{ $pop }->{ $chrom }->{ 'n' };
	
		print $out "$pop\t$chrom\t$n";
	
		foreach my $other_pop ( sort keys %{ $min_fst_table } ) {
		
			my $n_min = $min_fst_table->{ $pop }->{ $chrom }->{ 'group' }->{ $other_pop };
			
			print $out "\t" , $n_min / $n;
		
		}
		
		say $out "";
	
	}

}
