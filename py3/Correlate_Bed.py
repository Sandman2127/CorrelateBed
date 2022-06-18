#!/usr/bin/env python

import argparse
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
plt.switch_backend('agg') # makes sure matplotlib doesn't try to plot in the foreground


parser = argparse.ArgumentParser(description='Correlate the region values in column 5 of a bed file against another bedfile. i.e. H3K36me3 vs DNA methylation ')
parser.add_argument('--bed1', type=str,default=False,help='A normal bedfile you want to correlate regions against bed2, requires columns: chromosome\tstart\tend\tanyValue\tscore', required=True)
parser.add_argument('--bed2', type=str,help='A normal bedfile you want to correlate regions against bed1, requires columns: chromosome\tstart\tend\tanyValue\tscore', required=True)

args = parser.parse_args()
bed1=args.bed1
bed2=args.bed2

def parse_beds():
    bedlist1 = []
    bedlist2 = []
    bed1count = 0
    bed2count = 0
    matches_conditions = 0
    with open(bed1,'r') as inF1:
        lines = inF1.readlines()
        for line in lines:
            line = line.split()
            chrom1, start1, end1, score1 = str(line[0]),int(line[1]),int(line[2]),float(line[4])
            bed1count += 1
            with open(bed2,'r') as inF2:
                linez = inF2.readlines()
                for lin in linez:
                    lin = lin.split()
                    chrom2, start2, end2, score2 = str(lin[0]),int(lin[1]),int(lin[2]),float(lin[4])
                    bed2count += 1
                    #assumes a perfect match between file 1 & 2 for chromosome,position and region
                    if chrom1 == chrom2 and start1 == start2 and end1 == end2:
                        #the append is an ordered append, adding the scores inline with one another for downstream correlations
                        bedlist1.append(score1)
                        bedlist2.append(score2)
                        matches_conditions += 1 

    combined_list = [bedlist1,bedlist2]
    print("[STDOUT]: Example Bedlists:")
    print("[STDOUT]: BedList1:",bedlist1[0:5],"\n[STDOUT]: Bedlist2:",bedlist2[0:5])
    print("[STDOUT]: Out of:",bed1count,"regions in --bed1",float(100*(matches_conditions/bed1count)),"% or",matches_conditions,"exact matching regions were found in --bed2")
    print("[STDOUT]: The correlation of the measured values in these regions are as follows:")
    return combined_list

def calculate_correlation_of_matching_regions(bedlist1,bedlist2):
    # calculate spearman's correlation
    corrp, _ = pearsonr(bedlist1, bedlist2)
    print('[STDOUT]: Pearsons correlation: %.5f' % corrp)
    # calculate spearman's correlation
    corrs, _ = spearmanr(bedlist1, bedlist2)
    print('[STDOUT]: Spearmans correlation: %.5f' % corrs)
    
    return(corrp)

def plot_correlation(corrp,blist1,blist2):
    plotName = bed1[:-4] + "_vs_" + bed2[:-4] +  " Pearsons corr:" + str(round(corrp,3))
    plt.scatter(blist1,blist2)
    plt.title(plotName,)
    plt.xlabel('--bed1 values')
    plt.ylabel('--bed2 values')
    plt.grid(False)
    plot_name_jpg = bed1[:-4] + "_vs_" + bed2[:-4] + "_correlation_scatter_plot.jpg"
    plot_name_eps = bed1[:-4] + "_vs_" + bed2[:-4] + "_correlation_scatter_plot.eps"
    plt.savefig(plot_name_jpg,dpi=300)
    plt.savefig(plot_name_eps,dpi=300)
    plt.clf()


def main():
    combined_list = parse_beds()
    bedlist1, bedlist2 = combined_list[0],combined_list[1]
    corrp = calculate_correlation_of_matching_regions(bedlist1,bedlist2)
    plot_correlation(corrp,bedlist1,bedlist2)

if __name__ == "__main__":
    main()


