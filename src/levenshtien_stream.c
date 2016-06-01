#include <string.h>
#include <stdio.h>

static int distance (const char * word1,
		int len1,
		const char * word2,
		int len2)
{
	int matrix[len1 + 1][len2 + 1];
	int i;
	for (i = 0; i <= len1; i++) {
		matrix[i][0] = i;
	}
	for (i = 0; i <= len2; i++) {
		matrix[0][i] = i;
	}
	for (i = 1; i <= len1; i++) {
		int j;
		char c1;

		c1 = word1[i-1];
		for (j = 1; j <= len2; j++) {
			char c2;

			c2 = word2[j-1];
			if (c1 == c2) {
				matrix[i][j] = matrix[i-1][j-1];
			}
			else {
				int delete;
				int insert;
				int substitute;
				int minimum;

				delete = matrix[i-1][j] + 1;
				insert = matrix[i][j-1] + 1;
				substitute = matrix[i-1][j-1] + 1;
				minimum = delete;
				if (insert < minimum) {
					minimum = insert;
				}
				if (substitute < minimum) {
					minimum = substitute;
				}
				matrix[i][j] = minimum;
			}
		}
	}
	return matrix[len1][len2];
}

int main ()
{
	const char * word1;
	const char * word2;
	int len1;
	int len2;
	int d;

	/* Read in barcodes from STDIN  
	This is the maximum number of barcodes to compare
	in a group of barcodes related by a prefix or suffix */
	char uids[100000][13];
	
	char uid_in[13];
	int l;  // To count the number of lines read and store into the UIDs array
	while ( fgets(uid_in, 13, stdin) != NULL ) {
		// Chomp
		size_t ln = strlen(uid_in) - 1;
		if (uid_in[ln] == '\n') uid_in[ln] = '\0';
		
		strcpy(uids[l], uid_in);
		l++;
		// printf("read number: %s\n", uid_in);
	}
		
	// word1 = "pieces";
	// word2 = "paces";
	
	int i;
	int j;
	for ( i = 0; i < (l-1); i++ ) {
		for ( j = (i+1); j < l; j++ ) {
	
			word1 = uids[i];
			word2 = uids[j];
		
			len1 = strlen (word1);
			len2 = strlen (word2);
			d = distance (word1, len1, word2, len2);
			if ( d == 1 ) {
				printf("EditDist%d:\t%s\t%s\n", d, word1, word2);
				fflush(stdout);
				}
		}
	}
				
	// printf ("The edit distance between %s and %s is %d.\n",
			// word1, word2, d);
	return 0;
}
