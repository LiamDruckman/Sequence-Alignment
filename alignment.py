#!/usr/bin/python3

import sys
import argparse

def main():
    # sequence1 = "TCTGGCGATACTT"
    # sequence2 = "AGACCGCTATGAA"
    args = process_command_line()

    # get the sequences from the command line
    sequence1 = args.sequence1
    sequence2 = args.sequence2

    # set gap penalty
    gap_penalty = args.gap

    # set whether to use a blosum matrix or not and the match and mismatch scores
    use_blosum = args.blosum
    match_score = args.match
    mismatch_score = args.mismatch

    # set whether to use local alignment or not
    local = args.local

    # set the other optional arguments
    show_matrix = args.matrix
    show_arrows = args.no_arrows
    show_alignment = args.no_alignment
    show_score = args.no_score

    # calculate the score matrix and predecessor matrix
    matrix, predecessor = calculate(sequence1, sequence2, gap_penalty, match_score, mismatch_score, local, use_blosum)

    # print the matrix
    if show_matrix:
        print_matrix(matrix, predecessor, sequence1, sequence2, show_arrows)

    alignments = find_alignment(matrix, predecessor, sequence1, sequence2, local)

    if show_alignment:
        # print all alignments
        for i, alignment in enumerate(alignments):
            alignment1, alignment_match, alignment2 = alignment
            spacing = "\t" * 3
            print(f"Alignment {i+1}:")
            print(spacing + alignment1)
            print(spacing + alignment_match)
            print(spacing + alignment2)
            print()

            if show_score:
                # print the score and percentage identity
                print_score(matrix, alignment_match, local)


# Function to perform the calculation. Returns the score matrix and the predecessor matrix
def calculate(sequence1, sequence2, gap_penalty, match_score, mismatch_score, local, use_blosum):
    # import the blosum matrix if necessary
    if use_blosum:
        import blosum as bl
        scoring_matrix = bl.BLOSUM(62)

    # create a matrix of zeros
    matrix = [[0.0 for _ in range(len(sequence1) + 1)] for _ in range(len(sequence2) + 1)]

    # create a predecessor matrix to keep track of where the score came from
    predecessor = [[[] for _ in range(len(sequence1) + 1)] for _ in range(len(sequence2) + 1)]

    # fill in the first row and column with gap penalties
    for i in range(len(sequence1) + 1):
        matrix[0][i] = float(gap_penalty * i) if not local else 0.0
        if not local: predecessor[0][i].append("insert")
    for j in range(len(sequence2) + 1):
        matrix[j][0] = float(gap_penalty * j) if not local else 0.0
        if not local: predecessor[j][0].append("delete")
    predecessor[0][0].clear()

    # fill in the rest of the matrix saving where the score could have come from
    for i in range(1, len(sequence1) + 1):
        for j in range(1, len(sequence2) + 1):

            # if using a blosum matrix, get the score from the matrix
            if use_blosum:
                match = matrix[j - 1][i - 1] + scoring_matrix[sequence1[i - 1]][sequence2[j - 1]]
                if match == float("-inf"):
                    sys.stderr.write("Error: Invalid amino acid sequence. Please check your sequences and try again.\n")
                    sys.exit(1)

            elif sequence1[i - 1] == sequence2[j - 1]:
                match = matrix[j - 1][i - 1] + match_score
            else:
                match = matrix[j - 1][i - 1] + mismatch_score
            delete = matrix[j - 1][i] + gap_penalty
            insert = matrix[j][i - 1] + gap_penalty

            score_list = [match, delete, insert]
            # Can't go negative in local alignment
            if local: score_list.append(0.0)

            matrix[j][i] = max(score_list)

            # save where the score came from
            if matrix[j][i] == match:
                predecessor[j][i].append("match")
            if matrix[j][i] == delete:
                predecessor[j][i].append("delete")
            if matrix[j][i] == insert:
                predecessor[j][i].append("insert")

    return matrix, predecessor


# Function to print the matrix in a nice format
def print_matrix(matrix, predecessor, sequence1, sequence2, show_arrows):
    # Set up the arrows
    up_arrow = "↑" if show_arrows else " "
    left_arrow = "←" if show_arrows else " "
    diagonal_arrow = "↖" if show_arrows else " "

    # Look at every value in the matrix and find the longest one using map and zip
    padding = max(map(len, [str(item) for row in matrix for item in row])) + 1

    # Print the matrix
    # Top row comes first
    print("     | " + f"{'∅':^{padding * 2}}", end="")
    for letter in sequence1:
        print(f"{letter:^{padding * 2}}", end="")
    print()
    print("_____|" + "_" * (padding * 2) * (len(sequence1) + 1))

    # Then the rest of the rows
    for i in range(len(sequence2) + 1):
        if i == 0:
            curr_line = "  ∅  |"
        else:
            curr_line = "  " + sequence2[i - 1] + "  |"

        # print the matrix with each box padded for the longest value
        prev_line = "     |"
        for j in range(len(sequence1) + 1):
            prev_line += (diagonal_arrow if "match" in predecessor[i][j] else " ") + " " * (padding - 1) + (up_arrow if "delete" in predecessor[i][j] else " ") + " " * (padding - 1)
            curr_line += (left_arrow if "insert" in predecessor[i][j] else " ") + f"{matrix[i][j]:^{padding * 2 - 1}}"

        print(prev_line)
        print(curr_line)

    print()


# Function to find the alignment. Returns the alignment
def find_alignment(matrix, predecessor, sequence1, sequence2, local):
    alignments = []
    
    def dfs(i, j, alignment1, alignment_match, alignment2):
        nonlocal alignments

        if (i == 0 and j == 0) or (local and matrix[j][i] == 0.0):
            alignments.append((alignment1[::-1], alignment_match[::-1], alignment2[::-1]))
            return

        if "match" in predecessor[j][i]:
            dfs(i - 1, j - 1, alignment1 + sequence1[i - 1], alignment_match + ("|" if sequence1[i - 1] == sequence2[j - 1] else " "), alignment2 + sequence2[j - 1])
        if "delete" in predecessor[j][i]:
            dfs(i, j - 1, alignment1 + "-", alignment_match + " ", alignment2 + sequence2[j - 1])
        if "insert" in predecessor[j][i]:
            dfs(i - 1, j, alignment1 + sequence1[i - 1], alignment_match + " ", alignment2 + "-")

    if local:
        top_score = max(map(max, matrix))
        coordinates = [(i, j) for i in range(len(matrix)) for j in range(len(matrix[i])) if matrix[i][j] == top_score]
        for (i, j) in coordinates:
            dfs(i, j, "", "", "")
    else:
        dfs(len(sequence1), len(sequence2), "", "", "")

    return alignments


# Function to calculate and print the score and percentage identity of the alignment
def print_score(matrix, alignment_match, local):
    score = max(map(max, matrix)) if local else matrix[-1][-1]
    print("Score:" + "\t" * 3 + str(score))
    print("Percent Identity:" + "\t" + str(round(alignment_match.count("|") / len(alignment_match) * 100, 2)) + "%")
    print()


# Function to process command line arguments
def process_command_line():
    parser = argparse.ArgumentParser(description='Alignment Program. Only uses a linear gap penalty. Uses the Needleman-Wunsch algorithm for global alignment and the Smith-Waterman algorithm for local alignment.')
    parser.add_argument('sequence1', metavar='sequence1', type=str, help='First sequence. Can be any sequence of characters unless a scoring matrix is chosen')
    parser.add_argument('sequence2', metavar='sequence2', type=str, help='Second sequence. Can be any sequence of characters unless a scoring matrix is chosen')
    parser.add_argument('--local', action='store_true', default=False, help='Choose to use local alignment. Global alignment is done by default')
    parser.add_argument('--blosum', action='store_true', default=False, help='Choose to use a Blosum62 scoring matrix. Default is to use a standard match or mismatch score. If chosen the provided sequences must be amino acid sequences')
    parser.add_argument('--matrix', action='store_true', default=False, help='Choose to show the matrix. Matrix is not shown by default')
    parser.add_argument('--no-arrows', action='store_false', default=True, help='Choose to have retracing arrows - only works if matrix is shown. Arrows are shown by default')
    parser.add_argument('--no-alignment', action='store_false', default=True, help='Choose to not show the alignment. Alignment is shown by default')
    parser.add_argument('--no-score', action='store_false', default=True, help='Choose to not show the score and percentage identity of the alignment - only works if alignment is shown. Score is shown by default')
    parser.add_argument('-g', '--gap', metavar='gap', type=int, help='Gap penalty. Default is -8', default=-8)
    parser.add_argument('-m', '--match', metavar='match', type=int, help='Match score. Only works if no scoring matrix is chosen. Default is 5', default=5)
    parser.add_argument('-s', '--mismatch', metavar='mismatch', type=int, help='Mismatch score. Only works if no scoring matrix is chosen. Default is -4', default=-4)
    args = parser.parse_args()
    return args


# Driver Code
if __name__ == '__main__':
    main()
