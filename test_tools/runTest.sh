

# This test script will create a ./test_results directory in 
# the Debarcer root directory and run the pipeline on 
# a small test file.
# 
# The test should complete in about 2 minutes.

# Change to the root directory if the test script is run
# while in the test_tools directory
if [ -e ../runDebarcer.sh ]
then
	cd ..
fi

if [ -e runDebarcer.sh ]
then

	mkdir -p test_results
	cd test_results

	ln -s ../test_tools/Sample_Test.R1.fastq.gz
	echo "Running Debarcer in "`pwd`
	
	runDebarcer.sh -r -f Sample_Test.R1.fastq.gz -n Sample_Test
	
fi	
