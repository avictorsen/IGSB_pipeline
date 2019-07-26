"""
    Need to install a more up-to-date varsion of bedtools before invoking Jupyter
      type: module load bedtools/2.21.0
	(1) filters bedgraph based on threshold;
	(2) merges adjacent basepairs that are over threshold;
        (3) retains peaks that satisfy min/max length criteria;
	(4) merges any peaks that are closer than the inter-peak distance cutoff -or-
            alternatively keeps just the highest peak (this is beta functionality)

    Max length is typically defaulted to be very large
    Outputs a bed file (default col4 is the sum of the bedgraph scores; sorted by chrom;start;stop)
    Generate ID: will auto generate a integer list as a ID number (1... number of peaks). This will
      be reported as column 4 and the bedgraph scores will be shifted to column 5 as per standard bed format
    Note the peak score for merged peak is the *just* the sum of the two individual peaks not the
      total score in the merged region (i.e. there could be some sub-threshold scores in the intervening
      space that won't be included)
    Assumes bedgraph in standard format <chr> <start> <stop> <score>
    Output_name = option for user defined name (type with '...'), otherwise will generate name bedgraph_peaks.bed
    Delete_overlap_bed = option to add path to bedfile (as string), whereby any peaks that overlap this bed file will be discarded

___________
Usage: py_peak_calling.py [-f <folder> -t <threshold> --min <min_length> -i <inter_peak_distance> -m <merge_close_peaks> -1 <keep_highest_close_peak> --max <max_length> --ids <generate_ID> -o <output_name> --del <delete_overlap_bed>]

Options:
    -h --help
    -f <file>                      Path to bedgraph file [default: None]
    -t <threshold>                 [default: 40.0]
    --min <min_length>             [default: 1]
    -i <inter_peak_distance>       [default: 1]
    -m <merge_close_peaks>         [default: False]
    -1 <keep_highest_close_peak>   [default: True]
    --max <max_length>             [default: 10000]
    --ids <generate_ID>            [default: True]
    -o <output_name>               [default: None]
    --del <delete_overlap_bed>     [default: None]

___________
"""

from docopt import docopt
import sys

args = docopt(__doc__)

# %load /home/pskene/bin/py_peak_calling.py
def py_peak_calling(bedgraph, threshold, min_length, inter_peak_distance, merge_close_peaks, keep_highest_close_peak, max_length, generate_ID, output_name, delete_overlap_bed):
    import pybedtools
    import glob
    from pybedtools import BedTool
    import pandas as pd
    import csv

    if merge_close_peaks==keep_highest_close_peak:
        print 'Exiting... merge_close_peaks and keep_highest_close_peak set the same'
        sys.exit()

    #generate name for output
    bedgraph_name = glob.glob(bedgraph)
    filtered_name = bedgraph_name[0].replace('.bedgraph', 'filtered.bedgraph')

    if output_name != 'None':
        filename = output_name

    elif output_name == 'None':
        filename = bedgraph_name[0].replace('.bedgraph', '_peaks.bed')

    print 'input bedgraph file: ' + bedgraph_name[0]
    print 'output filename: ' + filename

    #import data as BedTool
    data = BedTool(bedgraph)
    print 'total sites read: ',
    print len(data)

    #retains intervals above threshold
    above_thresh = data.filter(lambda b: float(b.name) >= float(threshold)).saveas(filtered_name)
    print 'sites above threshold: ',
    print len(above_thresh)
    if len(above_thresh) == 0:
        print 'no regions are above the threashold\n'
        sys.exit()

    #merge adjacent above threshold regions and sum bedgraph scores (assumes bedgraph score in col 4)
    #d max distance between merged peaks, c: column modified
    merge_regions=above_thresh.merge(d=10, c=4, o='sum').saveas('temp.bed')

    #filter based on length criteria
    peaks = BedTool(merge_regions.filter(lambda x: len(x) >= min_length and len(x) <= max_length)).saveas('temp2.bed')
    print 'number of regions identified: ' + str(peaks.count())

    if merge_close_peaks=='True':
        if len(peaks) == 0:
            #merge the bonafide peaks if they they are shorter than the inter peak distance and sum scores and sort
            print 'merging peaks that are closer than: ' + str(inter_peak_distance)
            merge_peaks = peaks.merge(d=inter_peak_distance, c= 4, o='sum').sort().saveas('temp3.bed')
        if len(peaks) > 0:
            print 'no regions can be merged'
            merge_close_peaks='False'
            keep_highest_close_peak='True'

    if keep_highest_close_peak=='True':
        #need to read each line to find close peaks and throw away the one with the lowest score out of the two
        print 'entering loop'
#        if len(peaks) > 0:
        peaks.saveas('temp_input.bed')

        #print 'before keeping highest, number of regions identified: ' + str(BedTool('temp_input.bed').count())

        last_line = [str(item) for item in (BedTool('temp_input.bed').to_dataframe().tail(n=1).iloc[0,:].tolist())]

        with open('temp_input.bed') as myfile:
            with open('test_output.bed', 'w') as output:
                file_output = csv.writer(output, delimiter='\t')

                prev_line = None

                for line in csv.reader(myfile, delimiter='\t'):
                     print 'testing line: ' +str(line)

                     if prev_line is None:
                         prev_line = line
                         print

                     elif float(prev_line[2])+float(inter_peak_distance) <= float(line[1]):
                         print 'prev_line: ' + str(prev_line)
                         print 'line: ' + str(line)
                         print 'features far apart, so adding'
                         print
                         file_output.writerow(prev_line)
                         prev_line = line

                     else:
                         print 'prev_line: ' + str(prev_line)
                         print 'line: ' + str(line)
                         print 'features must be close'
                         print
                     if float(prev_line[3]) < float(line[3]):
                         prev_line = line
                         print 'prev_line smaller, so new prev_line'
                         print 'prev_line: ' + str(prev_line)
                         print

                print 'finished reading lines'
                print line
                print last_line
                if line==last_line:
                    print 'must be last line'
                    file_output.writerow(prev_line)

        merge_peaks = BedTool('test_output.bed')
    sys.exit()
    print 'number of peaks found: ' + str(merge_peaks.count())

    if delete_overlap_bed!=None:
        print 'delete_overlap_bed provided: ' + delete_overlap_bed
        merge_peaks = merge_peaks.intersect(b=delete_overlap_bed, v=True)
        print 'number of peaks retained: ' + str(merge_peaks.count())

    if not generate_ID:
        print 'saving sorted peak bed file with no ID'

        merge_peaks.saveas(filename)

    if generate_ID:
        print 'saving sorted peak bed file with ID names'

        #change to pandas dataframe
        DF_peaks = merge_peaks.to_dataframe()

        #insert new column with id: 1.... # of peaks
        DF_peaks.insert(3, 'id', ['id' + str(item) for item in range(1, (len(DF_peaks)+1))])

        ['id' + str(item)  for item in range(1, 5)]
        #save output
        DF_peaks.to_csv(filename, sep = '\t', header = False, index = False)

    return 'Finished'


py_peak_calling(bedgraph=args['-f'], threshold=args['-t'], min_length=args['--min'], inter_peak_distance=args['-i'], merge_close_peaks=args['-m'], keep_highest_close_peak=args['-1'], max_length=args['--max'], generate_ID=args['--ids'], output_name=args['-o'], delete_overlap_bed=args['--del'])
