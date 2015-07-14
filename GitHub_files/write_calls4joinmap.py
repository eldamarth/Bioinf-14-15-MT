#!usr/bin/env python
## Generate call files for mapping in Joinmap
## Input is SNPCx.pl output with all clustered and core SNPs
## Output is call files suitable to JoinMap

print('USAGE: write start(basepath, [lgnum], [perc]) to execute script.')

def body():
    
    with open("clusterout.txt") as f1:
        f = f1.readlines()

    coresnp = set()    # set([coresnps])
    assochgs = {}   # dict[snp][HG] = 1
    associnf = {}   # dict[snp][assoc] = phase
    hgs = {}
    for count, line in enumerate(f, 1):
        try:
            if not '/' in line:
                coresnp.add(line.split()[0])
                hgs[line.split()[1]] = {}
            else:
                lineinf = line.rstrip().split()
                if not lineinf[2] in assochgs:
                    assochgs[lineinf[2]] = {}
                    associnf[lineinf[2]] = {}
                if lineinf[1] == 'NORMAL':
                    lineinf[1] = 1
                elif lineinf[1] == 'INVERS':
                    lineinf[1] = -1
                assochgs[lineinf[2]][lineinf[0].split('/')[1]] = 1
                associnf[lineinf[2]][lineinf[0].split('/')[0]] = lineinf[1]
        except:
            with open('errorlog.txt', 'w') as err:
                err.writelines(str(count) + '\t' + line)

    difhg = set()      # Var to store snps mapped to >1 HG > difHG.txt
    for snp in assochgs:
        if len(assochgs[snp]) > 1:
            difhg.add(snp)


    def associatesnps(snp):
        for snp2 in associnf[snp]:
            if not snp2 in phase:
                phase[snp2] = phase[snp] * associnf[snp][snp2]
                if snp2 in associnf:
                    associatesnps(snp2)
            
    phase = {}  # Phase of core snps regarding their associations
    for snp in coresnp:
        if not snp in phase:
            if not snp in associnf:
                phase[snp] = 1   # NORMAL by default
            else:
                phase[snp] = 1   # Seed snp set to NORMAL
                associatesnps(snp)  # Recursive search and phasing of core snps
    

    snps = []
    classif = '(a,h,b)'
    phasing = []
    outsnp = difhg.copy().union(coresnp)
    for line in f:
        snpinf = line.rstrip().split()
        if '/' in snpinf[0] and snpinf[2] not in outsnp:  #Associated SNPs:assoc/HG,phase,id,calls
            id_ = snpinf[2]
            core, homg = snpinf[0].split("/")
            if core in difhg:
                continue
            calls = [call.replace('-1','-') for call in snpinf[3:]]
            if phase[core] * associnf[id_][core] == 1:  #NORM*NORM,INV*INV
                genot = [x.translate({ord('0'): 'a', ord('1'): 'h',
                                      ord('2'): 'b'}) for x in calls]
                calline = '\t'.join(['A_' + id_ + '_' + homg.replace('-','')
                                     ,classif] + genot)+'\n'
                phasing.append(id_ + '\tNORMAL\n')
            elif phase[core] * associnf[id_][core] == -1:
                genot = [x.translate({ord('0'): 'b', ord('1'): 'h',
                                      ord('2'): 'a'}) for x in calls]
                calline = '\t'.join(['A_' + id_ + '_' + homg.replace('-','')
                                     ,classif] + genot)+'\n'
                phasing.append(id_ + '\tINVERS\n')
            snps.append(calline)
            outsnp.add(id_)
            hgs[homg]['assoc'] = hgs[homg].get('assoc',0) + 1

        elif not '/' in snpinf[0]:  #Core SNPs: id,pos,hg,calls
            calls = [call.replace('-1','-') for call in snpinf[3:]]
            if phase[snpinf[0]] == 1:
                genot = [x.translate({ord('0'): 'a', ord('1'): 'h',
                                      ord('2'): 'b'}) for x in calls]
                phasing.append(snpinf[0] + '\tNORMAL\n')
            elif phase[snpinf[0]] == -1:
                genot = [x.translate({ord('0'): 'b', ord('1'): 'h',
                                      ord('2'): 'a'}) for x in calls]
                phasing.append(snpinf[0] + '\tINVERS\n')
            name = str(round(eval(snpinf[2])))+'_' + snpinf[0] + snpinf[1].replace('-','')
            snps.append('\t'.join([name, classif] + genot)+'\n')
            hgs[snpinf[1]]['core'] = hgs[snpinf[1]].get('core',0) + 1

    return difhg, hgs, snps, coresnp, phasing

def start(basepath, lgnum='7', perc='90'):  # basepath/LGlgnum/perc
    import os
    import os.path

    for lgs in range(1, int(lgnum) + 1):
        lg = 'LG' + str(lgs)
        path = os.path.join(basepath, lg, perc)
        os.chdir(path)
        difhg, hgs, snps, coresnp, phasing = body()
        lgfile = lg + perc + '.txt'
        with open("unclusterout.txt") as f1:
            header = f1.readline()
            header = header.split()
            header.insert(1,'Classif')
            header = '\t'.join(header) + '\n'
            unclus = len(f1.readlines())
    
        with open(lgfile, 'w') as out:
            out.writelines(header)
            out.writelines(snp for snp in snps)

        with open("difHG.txt", 'w') as out2:
            out2.writelines('probeset_id\n')
            out2.writelines((snp + '\n') for snp in difhg)

        with open("stats.txt", 'w') as out3:
            out3.writelines(lg + ', perc ' + perc + '\n')
            out3.writelines('Total SNPs:\t' + str(len(snps) + len(difhg)) + '\n')
            out3.writelines('Clustered Core SNPs:\t' + str(len(coresnp))
                            + '\t(' + str(round(len(coresnp)/len(snps)*100,1)) + '%)\n')
            out3.writelines('Unclustered Core SNPs:\t' + str(unclus) + '\n')
            out3.writelines('SNPs mapping in different homeologs:\t' + str(len(difhg))
                            + '\t(' + str(round(len(difhg)/len(snps)*100,1)) + '%)\n')
            out3.writelines('Remaining SNPs:\t'+ str(len(snps)) + '\n')
            out3.writelines('\n\n')
            out3.writelines('HG\tCore\tAssoc\tTotal\t%_LG\n')

            hgcore = 0
            hgassoc = 0
            hgtotal = 0
            for hg in sorted(hgs.keys()):
                hgcore = hgcore + hgs[hg]['core']
                hgassoc = hgassoc + hgs[hg]['assoc']
                hgtotal = hgtotal + round((hgs[hg]['core']+hgs[hg]['assoc'])/len(snps)*100,1)
                out3.writelines('\t'.join([hg, str(hgs[hg]['core']), str(hgs[hg]['assoc']),
                                str(hgs[hg]['core']+hgs[hg]['assoc']),
                                str(round((hgs[hg]['core']+hgs[hg]['assoc'])/len(snps)*100,1))])
                                + '\n')
            out3.writelines('\t'.join(['Total', str(hgcore), str(hgassoc), str(hgcore + hgassoc),
                            str(round(hgtotal))]) + '\n')

        with open("phasing.txt", 'w') as out4:
            out4.writelines('probeset_id\tphase\n')
            out4.writelines(snp for snp in phasing)
