#!/usr/bin/env perl

use strict;
use warnings;

# Alignment of two sequences by recursively calling the align function
# This is very inefficient, but should be quite easy to understand
# Huan Truong

my $match = +4;
my $mismatch = -5;
my $gap = -7;


# Align two sequeences
# This function takes in four parameters
# align($seq_x, $seq_y, $n_x, $n_y)
# seq_x and seq_y are the two sequences you want to get aligned
# n_x and n_y are the character index of the two sequences you want to get aligned
# it will return two things (score, alignment_x, alignment_y)

sub align {
	my ($seq_x, $seq_y, $n_x, $n_y) = @_;

	if ($n_x == -1) {
		return ($gap * ($n_y + 1), "-" x ($n_y + 1), substr($seq_y, 0, $n_y + 1));

	} elsif ($n_y == -1) {
		return ($gap * ($n_x + 1), "-" x ($n_x + 1), substr($seq_x, 0, $n_x + 1));

	} else {

		# we need to account for whether this is a match of mismatch
		my $base_on_left = substr($seq_x, $n_x, 1);
		my $base_on_top = substr($seq_y, $n_y, 1);

		my ($score_left, $align_left_x, $align_left_y) = align($seq_x, $seq_y, $n_x - 1, $n_y);

		$score_left = $score_left + $gap;

		my ($score_top, $align_top_x, $align_top_y) = align($seq_x, $seq_y, $n_x, $n_y - 1);

		$score_top = $score_top + $gap;


		my ($score_diag, $align_diag_x, $align_diag_y) = align($seq_x, $seq_y, $n_x - 1, $n_y - 1);

		if ($base_on_left eq $base_on_top) {
			$score_diag = $score_diag + $match;
		} else {
			$score_diag = $score_diag + $mismatch;
		}

		my $score_max = $score_diag;
		my $align_x = $align_diag_x . $base_on_left;
		my $align_y = $align_diag_y . $base_on_top;


		# The part over here selects the best score among the score_diag, left and top
		if ($score_top > $score_max) {
			$score_max = $score_top;
			$align_x = $align_top_x . "-";
			$align_y = $align_top_y . $base_on_top;
		}

		if ($score_left > $score_max) {
			$score_max = $score_left;
			$align_x = $align_left_x . $base_on_left;
			$align_y = $align_left_y . "-";
		}

		return ($score_max, $align_x, $align_y);
	}
}

# Align two sequeences
# This function takes in four parameters
# align($seq_x, $seq_y, $n_x, $n_y)
# seq_x and seq_y are the two sequences you want to get aligned
# n_x and n_y are the character index of the two sequences you want to get aligned
# It returns the best score possible between the two sequences

sub align_score {
	my ($seq_x, $seq_y, $n_x, $n_y) = @_;

	if ($n_x == -1) {
		return $gap * ($n_y + 1);

	} elsif ($n_y == -1) {
		return $gap * ($n_x + 1);

	} else {

		# we need to account for whether this is a match of mismatch
		my $base_on_left = substr($seq_x, $n_x, 1);
		my $base_on_top = substr($seq_y, $n_y, 1);

		my $score_left = align_score ($seq_x, $seq_y, $n_x - 1, $n_y);

		$score_left = $score_left + $gap;

		my $score_top = align_score ($seq_x, $seq_y, $n_x, $n_y - 1);

		$score_top = $score_top + $gap;

		my $score_diag = align_score ($seq_x, $seq_y, $n_x - 1, $n_y - 1);

		if ($base_on_left eq $base_on_top) {
			$score_diag = $score_diag + $match;
		} else {
			$score_diag = $score_diag + $mismatch;
		}

		# The part over here selects the best score among the score_diag, left and top
		my $score_max = $score_diag;

		if ($score_top > $score_max) {
			$score_max = $score_top;
		}

		if ($score_left > $score_max) {
			$score_max = $score_left;
		}

		return $score_max;
	}
}


sub align_score_fast {
	my ($seq_x, $seq_y, $n_x, $n_y) = @_;

	my @matrix;

	for (my $i = 0; $i < $n_x + 1; $i++) {
		$matrix[$i][0] = $gap * $i;
	}

	for(my $j = 0; $j < $n_y + 1; $j++) {
		$matrix[0][$j] = $gap * $j;
	}

	for (my $x = 1; $x < $n_x + 1; $x++) {
		for (my $y = 1; $y < $n_y + 1; $y++) {
			my $base_on_left = substr($seq_x, $x - 1, 1);
			my $base_on_top = substr($seq_y, $y - 1, 1);


			my $score_left = $matrix [$x - 1][$y];
			$score_left = $score_left + $gap;

			my $score_top = $matrix [$x][$y - 1];
			$score_top = $score_top + $gap;

			my $score_diag = $matrix [$x - 1][$y - 1];

			if ($base_on_left eq $base_on_top) {
				$score_diag = $score_diag + $match;
			} else {
				$score_diag = $score_diag + $mismatch;
			}

			my $score_max = $score_diag;

			if ($score_top > $score_max) {
				$score_max = $score_top;
			}

			if ($score_left > $score_max) {
				$score_max = $score_left;
			}

			$matrix[$x][$y] = $score_max;
			#print ("$x $y $score_max\n")
		}
	}

	return $matrix[$n_x][$n_y];
}



my $seq_1 = "TGAGCTAGTA";
my $seq_2 = "GAGCTAAA";

#my $alignment_score = align_score($seq_1, $seq_2, length($seq_1) - 1, length($seq_2) - 1);
#print ("Score of the alignment of $seq_1 and $seq_2 is $alignment_score\n");

my $alignment_score_fast = align_score_fast($seq_1, $seq_2, length($seq_1), length($seq_2));
print ("Score of the alignment of $seq_1 and $seq_2 is $alignment_score_fast\n");

my ($alignment_score_full, $alignment_x, $alignment_y) = align($seq_1, $seq_2, length($seq_1) - 1, length($seq_2) - 1);
print ("Score of the alignment of $seq_1 and $seq_2 is $alignment_score_full\n$alignment_x\n$alignment_y\n" );

