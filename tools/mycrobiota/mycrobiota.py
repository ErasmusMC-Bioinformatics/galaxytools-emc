import os
import argparse
import csv
import math
from subprocess import call


def main():
    print("Welcome to the MYcrobiota suite")
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--command', required=True,
                        help="What action to perform")
    parser.add_argument('-ct', '--count_table', action='append',
                        help="mothur count table")
    parser.add_argument('-cp', '--copies', help="copies of NC for samples")
    parser.add_argument('-nccp', '--nc_copies', help="copies of NC for itself")
    parser.add_argument('-t', '--taxonomy', action='append',
                        help="mothur taxonomy file")
    parser.add_argument('-s', '--shared_file', action='append',
                        help="mothur shared file")
    parser.add_argument('-otu', '--otutable', action='append',
                        help="mothur OTU table")
    parser.add_argument('-f', '--fasta', action='append', help="fasta")
    parser.add_argument('-sl', '--summary_log', action='append',
                        help="mothur summary log file")
    parser.add_argument('-o', '--outfile', help="output file")
    parser.add_argument('-od', '--outdir', help="output directory", default="")
    parser.add_argument('-lv', '--level', help="taxonomy level")
    parser.add_argument('-nc', '--negative_control',
                        help="sample name of the negative control")
    parser.add_argument('-ncs', '--negative_control_species',
                        help="species name of the negative control",
                        default="Oscillatoria")
    parser.add_argument('-r', '--replicate_suffix',
                        help="suffix to identify replicates")
    parser.add_argument('-l', '--label', action='append',
                        help="label for count table")
    parser.add_argument('--with-otu', dest='with_otu', action='store_true',
                        default=False)
    args = parser.parse_args()

    try:
        os.mkdir(args.outdir)
    except OSError:
        pass

    print("Running command: "+args.command)

    if args.command == 'counttable_totals':
        count_table_totals(args.count_table[0], args.outdir, args.outfile)

    elif args.command == 'qc_report':
        if args.count_table:
            qc_report(args.count_table, args.label, 'counttables', args.outdir)
        elif args.summary_log:
            qc_report(args.summary_log, args.label, 'summarylogs', args.outdir)

    elif args.command == 'create_krona_plot':
        create_krona_plot(args.taxonomy, args.outdir, args.with_otu)

    elif args.command == 'create_krona_plot_multisample':
        create_krona_plot_multisample(args.taxonomy, args.shared_file,
                                      args.level, args.outdir, args.with_otu)

    elif args.command == 'correct_replicates':
        correct_replicates(args.shared_file, args.taxonomy, args.outdir,
                           args.replicate_suffix, args.copies,
                           args.negative_control, args.nc_copies,
                           args.negative_control_species)

    elif args.command == 'make_multi_otutable':
        make_multi_otutable(args.taxonomy, args.shared_file, args.level,
                            args.outdir)

    elif args.command == 'otutable_add_blast_links':
        otutable_add_blast_links(args.otutable, args.fasta)

    elif args.command == 'split_multi_otutable':
        split_multi_otutable(args.otutable)

    else:
        print("unknown command. exiting")


def make_url(seq, baseurl):
    return baseurl+"?DATABASE=nr&PERC_IDENT=97&EXCLUDE_SEQ_UNCULT=on&" \
                   "HITLIST_SIZE=10&FILTER=L&FILTER=m&FILTER=R&EXPECT=10&" \
                   "FORMAT_TYPE=HTML&PROGRAM=blastn&CLIENT=web&" \
                   "SERVICE=megablast&PAGE=Nucleotides&CMD=Put&QUERY=" \
                  + seq.lower()


def make_RIDlink(RID, baseurl):
    return "<a target=\"_blank\" href=\""+baseurl+"?CMD=Get&RID="\
           + RID + "\">view results</a>"


def make_rerun_link(seq, baseurl):
    return "<a target=\"_blank\" href=\"" + baseurl +\
           "?DATABASE=nr&EXCLUDE_SEQ_UNCULT=yes&FILTER=L&FORMAT_TYPE=HTML" \
           "&PROGRAM=blastn&CLIENT=web&SERVICE=megablast&PAGE=Nucleotides&" \
           "CMD=Web&QUERY=" + seq.lower() + "\">send to BLAST</a>"


def otutable_add_blast_links(otutable, otureps):
    baseurl = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"

    # for each otu create blast search of corresponding representative sequence
    reps = [line for line in open(otureps[0], "r")]

    seqs = [r.rstrip('\n').replace('-', '') for r in reps if '>' not in r]
    seq_names = [r for r in reps if '>' in r]
    otulines = [line for line in open(otutable[0], "r")]

    # Add RID link and rerun link to table
    with open("otutable_with_blast.tsv", "w+") as outfile, open("filtered_otureps.fasta", "w+") as repsout:
        outfile.write(otulines[0].rstrip() + "\tBLAST\n")

        for otuline in otulines[1:]:
            otu = otuline.split('\t')[0]
            for i, seq in enumerate(seq_names):
                if otu in seq:
                    outfile.write(otuline.rstrip() + "\t" +
                                  make_rerun_link(seqs[i], baseurl) + "\n")
            # output otureps for these otus
            for i, seq in enumerate(reps):
                if otu in seq:
                    repsout.write(reps[i])
                    repsout.write(reps[i+1])


def summarylog_total(infile):
    with open(infile) as f:
        summarylog = f.readlines()
        for line in summarylog:
            if line.startswith('# of Seqs:') \
                    or line.startswith('total # of seqs:'):
                return int(line.split('\t')[-1])
    return None


def count_table_totals(infile, outdir='', outfile=''):
    """
    Given a Mothur counttable, calculate the total number of sequences for each
    sample. This can be appended as additional row in the count table by
    providing a file name.

    :param infile: Mothur count table.
    :param outfile: Optional. Write the count table with an additional row with
                    totals to this file
    :param outdir: Optional. Write output do this directory
    :return: A list with totals for all columns (samples) in the count table
    """
    # assume a single input file for now
    out_rows = []
    with open(infile) as f:
        count_table = csv.reader(f, delimiter='\t')

        header = next(count_table)
        totals = [0] * (len(header)-1)

        out_rows.append(header)
        for row in count_table:
            if row[0] != 'total':
                out_rows.append(row)
                for i in range(1, len(row)):
                    totals[i-1] += int(row[i])

        out_rows.append(["total"] + map(str, totals))

    # write count table with totals to file if requested
    if outfile:
        write_output(outdir, outfile, out_rows)

    return totals


def qc_report(infiles, label, inputtype, outdir=''):
    """
    Construct QC table from multiple count files
      - Report the number of sequences lost at each consecutive QC step
      - Create a multi-sample report and a separate report for each sample

    :param infiles: set of count tables
    :param label: labels for each step
    :param outdir: directory to place output files. Default: cwd
    :return:
    """
    assert len(infiles) == len(label), \
        "number of files and labels unequal, stopping"

    print("qcreport")
    previous_totals = []
    outlines = []
    lasttotal = None

    for (i, lab) in zip(infiles, label):
        with open(i, 'rb') as f:
            count_file = csv.reader(f, delimiter='\t')

            if inputtype == 'summarylogs':
                print("summarylogs")
                if not outlines:
                    outlines = [['Step', 'Total', 'Removed', 'Percentage']]

                # get total count
                total = summarylog_total(i)

                # calculate difference with last
                if not lasttotal:
                    outlines.append([lab, total, None, None])
                else:
                    diff = total - lasttotal
                    perc = float(diff)/float(lasttotal)*100.0
                    outlines.append([lab, total, diff, str("%.1f" % perc)+"%"])

                lasttotal = total

            else:
                # add header line to output
                if not outlines:
                    outlines = [['step'] + next(count_file)[1:]]

                # calculate totals of each column in the count file
                totals = count_table_totals(i)

                # calculate difference with previous count file
                if not previous_totals:
                    diffs = [""] * len(totals)
                else:
                    diffs = [" (" + str(t1-t2)+"; "
                             + str("%.1f" % (float(t1-t2)/float(t2)*100.0)) +
                             "%)" for t1, t2 in zip(totals, previous_totals)]

                outlines.append([lab] +
                                [str(a)+b for a, b in zip(totals, diffs)])
                previous_totals = totals

            # write multi-sample output file
            write_output(outdir, 'all_qctable.tsv', outlines)

            # write per-sample output files
            for j in range(2, len(outlines[0])):
                sample = outlines[0][j]
                sample_outlines = [[outlines_line[0], outlines_line[j]]
                                   for outlines_line in outlines]
                write_output(outdir, 'persample_qctable_'+sample+'.tsv',
                             sample_outlines)


def column(matrix, i):
    return [row[i] for row in matrix]


def mean(data):
    return sum(data) / float(len(data))


def stdev(data):
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    n = len(data)
    return math.sqrt(ss/(n-1))


def write_output(outdir, filename, outlines):
    with open(os.path.join(outdir, filename), 'wb') as of:
        out_table = csv.writer(of, delimiter='\t', lineterminator='\n')
        for row in outlines:
            out_table.writerow(row)


def correct_replicates(shared, taxonomy, outdir, replicate_suffix,
                       sample_copies, negative_control='', nc_copies=-1,
                       negative_control_species='Oscillatoria'):
    with open(shared[0], 'rb') as f, open(taxonomy[0], 'rb') as f2:
        shared_file = csv.reader(f, delimiter='\t')
        taxonomy_file = csv.reader(f2, delimiter='\t')

        # determine which OTU number is the control, Oscillatoria by default
        # (Bacteria;Cyanobacteria;Cyanobacteria;..;Oscillatoria)
        try:
            line = next(taxonomy_file)
            while negative_control_species not in line[2]:
                line = next(taxonomy_file)
            otu = line[0]
        except StopIteration:
            print("negative control species not found in taxonomy, Exiting")
            return 1


        ''' Calculate Copies '''
        # per replicate of sample and NC, determine correction factor,
        # (number Oscillatoria seqs/known copies of it)
        # correct all sequence counts with that
        myshared = [row for row in shared_file if row]
        newshared = [myshared[0]]
        newshared2 = [myshared[0]]
        newshared3 = [myshared[0]]
        oscil_column = myshared[0].index(otu)

        for row in myshared[1:]:
            if row[1].startswith(negative_control):
                copies = nc_copies
            else:
                copies = sample_copies

            correction_factor = float(row[oscil_column]) / float(copies)

            new_row = row[0:3]
            for count in row[3:]:
                new_row.append(float(count) / correction_factor)
            newshared.append(new_row)

        ''' Average copy counts across replicates  '''
        levels = set([row[0] for row in newshared[1:]])
        samples = set([row[1].split(replicate_suffix)[0]
                       for row in newshared[1:]])

        for level in levels:
            for sample in samples:
                neg = True if sample.startswith(negative_control) else False
                replicates = [row for row in newshared if row[0] == level
                              and row[1].startswith(sample)]
                num_otus = int(replicates[0][2])+3
                total = replicates[0][2]
                avg = [level, sample, total]

                for i in range(3, num_otus):
                    counts = column(replicates, i)
                    avg.append(mean(counts)) if 0 not in counts or neg \
                        else avg.append(0)

                newshared2.append(avg)

        ''' Correct for background '''
        # for each otu, subtract 3 times the standard deviation of
        # the negative control sample
        for level in levels:
            NC = [row for row in newshared if row[0] == level
                  and row[1].startswith(negative_control)]
            samples = [row for row in newshared2 if row[0] == level
                       and not row[1].startswith(negative_control)]
            num_otus = int(samples[0][2])+3

            for i in range(3, num_otus):
                m = mean(column(NC, i))
                sd = stdev(column(NC, i))
                corr = m + 3*sd

                for s in samples:
                    s[i] = max(0, int(round(s[i] - corr)))

            newshared3 += samples

        # remove Negative control species otu from table
        for i, row in enumerate(newshared3):
            del row[oscil_column]
            if i > 0:
                row[2] = int(row[2]) - 1

        # sort file or other mothur tools may segfault :/
        newshared3 = [newshared3[0]] + sorted(newshared3[1:],
                                              key=lambda a_entry: a_entry[0]
                                              if a_entry[0] != 'unique' else 0)

        f2.seek(0)
        taxonomy_out = [row for row in taxonomy_file if row and row[0] != otu]
        write_output(outdir, 'taxonomy_corrected.tsv', taxonomy_out)
        write_output(outdir, 'shared_corrected.tsv', newshared3)


def make_multi_otutable(taxonomy_file, shared_file, level, outdir):
    """
    Create an otu table from shared file and taxonomy file

    example output:

    OTU    sample1 sample2 .. sampleX Kingdom  Phylum  Class Order Family Genus
    Otu001 13      8       .. 91      Bacteria Bacteroidetes Bacteroidia ..
    ...

    :param taxonomy_file:
    :param shared_file:
    :param level:

    :return:
    """

    outlines = []
    samples = []
    # create multisample taxonomy from counts in shared file
    with open(taxonomy_file[0], 'r') as tax, open(shared_file[0]) as sh:
        taxonomy = csv.reader(tax, delimiter='\t')
        shared = csv.reader(sh, delimiter='\t')
        shared_header = next(shared)
        outlines.append(shared_header[3:])

        # get all taxonomies
        taxonomies = []
        for j, t in enumerate(taxonomy):
            if j > 0:
                taxonomies.append(filter(None, [x.split('(')[0]
                                                for x in t[2].split(';')]))

        for i, row in enumerate(shared):
            tax.seek(0)  # make sure to start at beginning of file each time
            if row[0] == level:
                samples.append(row[1])
                outlines.append(row[3:])

        transposed = map(list, zip(*outlines))
        header = ["OTU"] + samples + ["Kingdom", "Phylum", "Class", "Order",
                                      "Family", "Genus"]

        writelines = [header] + [a + b for a, b in zip(transposed, taxonomies)
                                 if a[1:] != ['0'] * len(a[1:])]

        # output corrected shared file
        write_output(outdir, "multi_otutable.tsv", writelines)


def split_multi_otutable(otutable, with_avg=True):
    fulltable = [line.strip().split('\t')
                 for line in open(otutable[0], 'r') if line]
    samples = [s.split('_')[0] for s in fulltable[0][1:-6]]
    numcols = len(fulltable[0])
    numreplicates = (numcols - 7) / len(set(samples))

    for sample in set(samples):
        outlines = []
        cols = [0] + [i+1 for i, s in enumerate(samples) if sample in s] \
            + [i for i in range(numcols-6, numcols)]
        for i, line in enumerate(fulltable):
            out = [line[j] for j in cols]
            if out[1:-6] != ['0'] * numreplicates:
                out.insert(-6, 'mean' if i == 0
                           else int(round(mean(map(int, out[1:-6])))))
                outlines.append(out)

        write_output('.', sample+'.otutable', outlines)


def create_krona_plot_multisample(taxonomy_file, shared_file, level, outdir,
                                  with_otu):
    """
    Create krona plots from a multisample taxonomy plot and a shared file.
    Create one multisample plot and a plot per individual sample

    :param taxonomy_file:
    :param shared_file:
    :param level: which level to use, e.g. unique/0.03/..
    :param with_otu:
    :return:
    """

    taxonomies = []

    # create taxonomy file per sample
    with open(taxonomy_file[0], 'r') as tax, open(shared_file[0]) as sh:
        taxonomy = csv.reader(tax, delimiter='\t')
        shared = csv.reader(sh, delimiter='\t')
        shared_header = next(shared)

        for i, row in enumerate(shared):
            tax.seek(0)  # make sure to start at beginning of file each time
            if row[0] == level:
                sample = row[1]

                outfile = os.path.join(outdir, sample+".tsv")
                taxonomies.append(outfile)
                with open(outfile, 'w+') as of:
                    out_table = csv.writer(of, delimiter='\t')
                    out_table.writerow(next(taxonomy))  # header line
                    for j, t in enumerate(taxonomy):
                        assert t[0] == shared_header[j+3], \
                            "OTU mismatch between taxonomy and shared file"
                        t[1] = row[j+3]
                        out_table.writerow(t + [shared_header[j+3]])

    # if one one sample in the shared file, don't create "allsamples" plot
    if len(taxonomies) == 2:
        taxonomies = taxonomies[1:]

    # make krona plot
    create_krona_plot(taxonomies, outdir, with_otu)


def create_krona_plot(taxonomy_files, outdir, with_otu):
    """
    Create a krona plot from one or more mothur taxonomy files

    :param taxonomy_files: mothur taxonomy file (output from classify.otu)
    :param outdir: directory to store krona-formatted outputs. Default=cwd
    :param with_otu: add OTU number as a level in the Krona plot? Default=True
    :return:
    """
    krona_input_files = []

    # convert taxonomy files to krona input.
    for tax in taxonomy_files:
        with open(tax, 'r') as f:
            taxonomy = csv.reader(f, delimiter='\t')
            out_rows = []

            next(taxonomy)  # skip header line
            for row in taxonomy:
                out_rows.append(
                    filter(None, [row[1]] + row[2].rstrip(";\n").split(';') +
                           [row[0] if with_otu else None]))

        outfile = os.path.join(outdir, tax.split("/")[-1]+"krona")
        krona_input_files.append(outfile)

        with open(outfile, 'w+') as f2:
            out_table = csv.writer(f2, delimiter='\t')
            for row in out_rows:
                out_table.writerow(row)

    # execute krona command
    call(["ktImportText"] + krona_input_files)


if __name__ == "__main__":
    main()
