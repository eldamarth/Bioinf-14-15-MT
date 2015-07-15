#!python3

import os
import os.path

def change_dir(path):
    import errno
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    os.chdir(path)

def file2set(file, field=1):  # 1-based field
    try:
        with open(file) as f:
            myset = set([line.rstrip().split()[field-1] for line in f])
        return myset
    except:
        print('Could not open file to obtain list')
        file = input('Please, reinsert file [path/]name: ')

def read_file(file, header=True, sep='\t'):
    try:
        comments = []
        with open(file) as f:
            for line in f:
                if line.startswith('#'):
                    comments.append(line)
                else:
                    break
            if header == True:
                header = line.rstrip().split(sep) 
                body = f.readlines()
                body = [line.rstrip().split(sep) for line in body]
                return comments, header, body
            else:
                body = [line.rstrip().split(sep)]
                body.extend([line.rstrip().split(sep) for line in f])
                return comments, body
    except:
        file = input('Could not open file, insert file [path/]name: ')
        read_file(file, header=header, sep=sep)

def write_output(joinedlist, outdir, outfilename):
    change_dir(outdir)
    with open(outfilename, 'w') as out:
        for line in joinedlist:
            out.writelines(line)

def single_dict(splitlist, key=0, store=1):
    _dict = {}
    for line in splitlist:
        _dict[line[key]] = line[store]
    return _dict

def write_snpfa(annotFile, outPath):
    import re
    fa = []
    com, header, annots = read_file(annotFile, sep=',')
    for line in annots:
        if 'snp' in line[9].lower() or 'codon' in line[9]:
            fa.append('>' + line[0].replace('"','') + '\n' + re.sub(
                     '\[.+?\]','N',line[6].replace('"','')) + '\n')
    write_output(fa, outPath, 'snp_seqs.fa')

def write_vcf(annotFile, outPath, outfile, sep='\t', flankpos=7):
    import re
    vcf = ['##fileformat=VCFv4.1\n']
    vcf.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                          'FILTER', 'INFO', 'FORMAT']) + '\n')
    com, header, ann = read_file(annotFile, sep=sep)
    for line in ann:
        var = re.split('[\[/\]]',line[flankpos])[1:3]
        vcf.append('\t'.join([line[3], line[4], line[0], var[0], 
                   var[1], '.', '.', '.']) + '\n')
    write_output(vcf, outPath, outfile)

def subfile(listfile, file, header=True, f1=1, f2=1, sep='\t', printcomments=False):
    row_names = file2set(listfile, field=f1)
    if header == True:
        comments, header, filebody = read_file(file, sep=sep)
        matched = [header] + [line for line in filebody if line[f2-1] in row_names]
    else:
        comments, filebody = read_file(file, sep=sep, header=False)
        matched = [line for line in filebody if line[f2-1] in row_names]
    if printcomments == True:
        return comments, matched
    return matched

def write_subfile(listfile, file, outPath, outFile, header=True, f1=1, f2=1, 
                  printcomments=False):
    if printcomments == True:
        comments, lines = subfile(listfile, file, header, f1, f2, 
                                  printcomments=printcomments)
        joinedlines = comments + ['\t'.join([l for l in line])+'\n' for line in lines]
    else:
        lines = subfile(listfile, file, header, f1, f2, printcomments=printcomments)
        joinedlines = ['\t'.join([l for l in line])+'\n' for line in lines]
    write_output(joinedlines, outPath, outFile)

def GCcalls2tassel(calls, header):
    tr = {'-1': '?:?', '0': 'A:A', '1': 'A:B', '2': 'B:B'}
    for line in calls:
        line[1:] = [tr[c] for c in line[1:]]
    tass_f = ['\t'.join(list(x)) + '\n' for x in zip(*[header]+calls)]
    return tass_f

def write_calls4tassel(callsfile, outpath=''):
    com, header, calls = read_file(callsfile)
    header = ['<marker>'] + header[1:]
    tass_f = GCcalls2tassel(calls, header)
    if not outpath:
        outpath = os.path.dirname(callsfile)
    write_output(tass_f, outpath, os.path.basename(callsfile).rsplit('.',1)[0]+'4Tass.txt')

def GCcalls2Joinmap(callsfile, outpath=''):
    header, calls = read_file(callsfile)[1:]
    header.insert(1, 'Classif')
    header = ['\t'.join([x for x in header]) + '\n']
    classif = '(a,h,b)'
    for line in calls:
        line[1:] = [call.replace('-1','-') for call in line[1:]]
        line[1:] = [x.translate({ord('0'): 'a', ord('1'): 'h',
                                 ord('2'): 'b'}) for x in line[1:]]
        line.insert(1, classif)
        header.append('\t'.join([x for x in line]) + '\n')
    if not outpath:
        outpath = os.path.dirname(callsfile)
    filename = os.path.basename(callsfile).rsplit('.',1)[0]+'4jm.txt'
    write_output(header, outpath, filename)

def write_polymorphic(callsfile, outfile, last_var=37):
    os.chdir(os.path.dirname(callsfile))
    polym = {'0', '1', '2'}
    with open(outfile, 'w') as out:
        with open(callsfile) as f:
            out.writelines(f.readline())
            for line in f:
                if polym.issubset(set(line.rstrip().split()[1:last_var])):
                    out.writelines(line)