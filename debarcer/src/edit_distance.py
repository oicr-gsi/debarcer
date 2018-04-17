
def edit_distance(a, b):
    """Returns the Hamming edit distance between two strings."""
    return sum(letter_a != letter_b for letter_a, letter_b in zip(a, b))