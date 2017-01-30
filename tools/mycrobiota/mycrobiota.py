import os
import argparse
import csv
from subprocess import call


def main():
    print("Welcome to MYcrobiota suite")
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--command', required=True, help="What action to perform")
    parser.add_argument('-ct', '--count_table', action='append', help="mothur count table")
    parser.add_argument('-t', '--taxonomy', action='append', help="mothur taxonomy file")
    parser.add_argument('-s', '--shared_file', action='append', help="mothur shared file")
    parser.add_argument('-sl', '--summary_log', action='append', help="mothur summary log file")
    parser.add_argument('-o', '--outfile', help="output file")
    parser.add_argument('-od', '--outdir', help="output directory", default="")
    parser.add_argument('-lv', '--level', help="taxonomy level")
    parser.add_argument('-nc', '--negative_control', help="sample name of the negative control")
    parser.add_argument('-r', '--replicate_suffix', help="suffix to identify replicates")
    parser.add_argument('-l', '--label', action='append', help="label for count table")
    parser.add_argument('--with-otu', dest='with_otu', action='store_true', default=False)
    args = parser.parse_args()

    try:
        os.mkdir(args.outdir)
    except OSError:
        pass

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
        create_krona_plot_multisample(args.taxonomy, args.shared_file, args.level, args.outdir, args.with_otu)

    elif args.command == 'correct_replicates':
        correct_replicates(args.shared_file, args.outdir, args.replicate_suffix, args.negative_control)

    elif args.command == 'shared_to_otutable':
        shared_to_otutable(args.taxonomy, args.shared_file, args.level, args.outdir)

    else:
        print("unknown command. exiting")


def summarylog_total(infile):
    with open(infile) as f:
        summarylog = f.readlines()
        for line in summarylog:
            if line.startswith('# of Seqs:') or line.startswith('total # of seqs:'):
                return int(line.split('\t')[-1])
    return None


def count_table_totals(infile, outdir='', outfile=''):
    """
    Given a Mothur count table, calculate the total number of sequences for each sample. This can be appended as
    additional     row in the count table by providing a file name.

    :param infile: Mothur count table.
    :param outfile: Optional. Write the count table with an additional row with totals to this file
    :param outdir: Optional. Write output do this directory
    :return: A list with the totals for all columns (samples) in the count table
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
    :param outdir: directory to place output files. Default: current working directory
    :return:
    """
    assert len(infiles) == len(label), "number of files and labels unequal, stopping"

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
                    diffs = [" ("+str(t1-t2)+"; "+str("%.1f" % (float(t1-t2)/float(t2)*100.0))+"%)"
                             for t1, t2 in zip(totals, previous_totals)]

                outlines.append([lab] + [str(a)+b for a, b in zip(totals, diffs)])
                previous_totals = totals

            # write multi-sample output file
            write_output(outdir, 'all_qctable.tsv', outlines)

            # write per-sample output files
            for j in range(2, len(outlines[0])):
                sample = outlines[0][j]
                sample_outlines = [[outlines_line[0], outlines_line[j]] for outlines_line in outlines]
                write_output(outdir, 'persample_qctable_'+sample+'.tsv', sample_outlines)


def column(matrix, i):
    return [row[i] for row in matrix]


def mean(data):
    return sum(data) / float(len(data))


def stdev(data):
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    n = len(data)
    return ss/(n-1)


def write_output(outdir, filename, outlines):
    with open(os.path.join(outdir, filename), 'wb') as of:
        out_table = csv.writer(of, delimiter='\t', lineterminator='\n')
        for row in outlines:
            out_table.writerow(row)


def correct_replicates(infile, outdir, replicate_suffix, negative_control=''):
    """
    Given a shared file, per sample, remove any OTUs not present in all replicates, and for those that remain, use
    average of counts of all replicates

    :param infile:
    :param outdir:
    :param replicate_suffix:
    :param negative_control:
    :return:
    """
    with open(infile[0], 'rb') as f:
        shared_file = csv.reader(f, delimiter='\t')
        out_lines = [next(shared_file)]  # header
        nc_stdevs = []

        # load all replicates of a sample (this assumes them to be in order in the file)
        peek = next(shared_file)
        while peek:
            replicates = []
            sample = peek[1].split(replicate_suffix)[0]
            while peek[1].split(replicate_suffix)[0] == sample:
                replicates.append(peek)
                try:
                    peek = next(shared_file)
                except StopIteration:
                    peek = False
                    break

            # calculate mean counts across replicates or set to zero if any of replicates had count zero
            averages = []
            stdevs = []
            if negative_control and replicates[0][1].startswith(negative_control):
                stdevs.append(replicates[0][0])  # add label to line

            for col in range(3, len(replicates[0])):
                counts = map(int, column(replicates, col))
                if 0 in counts:
                    averages.append(0)
                else:
                    averages.append(int(round(mean(counts))))

                # calculate normal control correction factor
                corr_factor = 0
                if negative_control and replicates[0][1].startswith(negative_control):
                    if 0 not in counts:
                        corr_factor = 3 * stdev(counts)
                    print replicates[0][0] + ': ' + replicates[0][1]
                    print counts
                    print corr_factor
                    stdevs.append(corr_factor)

            if stdevs:
                nc_stdevs.append(stdevs)

            # output single row per sample with corrected counts
            out_lines.append([replicates[0][0]] + [sample] + [replicates[0][2]] + averages)

        print nc_stdevs
        write_output(outdir, 'shared_dereplicated.tsv', out_lines)

        # if normal control sample was present, correct file with that
        if negative_control:
            # correct_negative_control('shared_dereplicated.tsv', outdir, negative_control)
            correct_negative_control('shared_dereplicated.tsv', outdir, negative_control, nc_stdevs)
        else:
            os.rename('shared_dereplicated.tsv', 'shared_corrected.tsv')


def correct_negative_control(infile, outdir, negative_control, stdevs):
    with open(os.path.join(outdir, infile)) as f:
        shared_file = csv.reader(f, delimiter='\t')
        corrected_lines = [next(shared_file)]  # header

        # per level in the shared file (unique, 0.03, ..), and per otu, calculate correction factor
        peek = next(shared_file)
        while peek:
            if peek[1] != negative_control:
                for level in stdevs:
                    if level[0] == peek[0]:
                        print peek
                        print level
                        # subtract correction factor from counts
                        corr_line = peek[0:3]
                        for i in range(3, len(peek)):
                            corr_line.append(int( max(0, int(peek[i]) - level[i-2])))

                        corrected_lines.append(corr_line)
                        print corr_line

            # get next line
            try:
                peek = next(shared_file)
            except StopIteration:
                peek = False
                break

        # output corrected shared file
        write_output(outdir, "shared_corrected.tsv", corrected_lines)


def correct_negative_control2(infile, outdir, negative_control):
    """
    If a negative control sample was added, correct the OTU counts in a shared file. Subtract 3 times the standard
    deviation of normal control sample counts from all counts in file.

    :param infile: shared file
    :param outdir: directory to output results
    :param negative_control: name of the normal control sample
    :return:
    """

    # calculate correction factor (3 * stdev(normal_control))
    with open(os.path.join(outdir, infile)) as f:
        shared_file = csv.reader(f, delimiter='\t')
        corrected_lines = [next(shared_file)]  # header
        correction = 1

        # per level in the shared file (unique, 0.03, ..), and per otu, calculate correction factor
        peek = next(shared_file)
        while peek:
            level_lines = []
            level = peek[0]
            while peek[0] == level:
                sample = peek[1]
                if sample == negative_control:
                    correction = 3 * stdev(map(int, peek[3:]))
                    print("correction factor: " + str(correction))
                else:
                    level_lines.append(peek)

                try:
                    peek = next(shared_file)
                except StopIteration:
                    peek = False
                    break

            # apply correction to samples
            for line in level_lines:
                corrected_lines.append(line[0:3] + map(lambda x: int(round(x-correction)), map(int, line[3:])))

    # output corrected shared file
    write_output(outdir, "shared_corrected.tsv", corrected_lines)


def shared_to_otutable(taxonomy_file, shared_file, level, outdir):
    """
    Create an otu table from shared file and taxonomy file

    example output:

    OTU    sample1 sample2 .. sampleX Kingdom  Phylum        Class       Order         Family         Genus
    Otu001 13      8       .. 91      Bacteria Bacteroidetes Bacteroidia Bacteroidales Prevotellaceae Prevotella
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

        # get all taxonomies
        taxonomies = []
        for j, t in enumerate(taxonomy):
            taxonomies.append(t[2].split(';'))

        for i, row in enumerate(shared):
            tax.seek(0)  # make sure to start at beginning of file each time
            if row[0] == level:
                samples.append(row[1])
                outlines.append(row[1:])

        transposed = map(list, zip(*outlines))

        header = ["OTU"] + samples + ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
        print header
        writelines = [header]
        writelines = [a+b for a, b in zip(transposed, taxonomies)]

        #print writelines
        # output corrected shared file
        write_output(outdir, "shared_with_taxonomy.tsv", writelines)


def create_krona_plot_multisample(taxonomy_file, shared_file, level, outdir, with_otu):
    """
    Create krona plots from a multisample taxonomy plot and a shared file. Create one multisample plot and a plot per
    individual sample

    :param taxonomy_file:
    :param shared_file:
    :param level: which level to use, e.g. unique/0.03/..
    :param with_otu:
    :return:
    """

    os.link(taxonomy_file[0], 'all-samples.tsv')
    taxonomies = ['all-samples.tsv']

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
                        assert t[0] == shared_header[j+3], "OTU mismatch between taxonomy and shared file"
                        t[1] = row[j+3]
                        out_table.writerow(t + [shared_header[j+3]])

    # make krona plot
    create_krona_plot(taxonomies, outdir, with_otu)


def create_krona_plot(taxonomy_files, outdir, with_otu):
    """
    Create a krona plot from one or more mothur taxonomy files

    :param taxonomy_files: mothur taxonomy file (output from classify.otu)
    :param outdir: Optional: directory to store krona-formatted outputs. Default=working directory
    :param with_otu: add OTU number as a level in the Krona plot? Default = True
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
                out_rows.append(filter(None, [row[1]] + row[2].rstrip(";\n").split(';') + [row[0] if with_otu else None]))

        outfile = os.path.join(outdir, tax.split("/")[-1])
        krona_input_files.append(outfile)

        with open(outfile, 'w+') as f2:
            out_table = csv.writer(f2, delimiter='\t')
            for row in out_rows:
                out_table.writerow(row)

    # execute krona command
    call(["ktImportText"] + krona_input_files)


if __name__ == "__main__":
    main()
