import pandas as pd

def InserInSilPCRlink(df,MaxL=600,genome='hg38'):
    if genome=='hg38':
        gnomadGenome = 'gnomad_r3'
    else:
        gnomadGenome = 'gnomad_r2_1'

    l_ucsc = f'https://genome.ucsc.edu/cgi-bin/hgPcr?org=Human&db={genome}&wp_target=genome&wp_f=___frwdP___&wp_r=___reverseP___&Submit=submit&wp_size={MaxL}&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0'
    
    lgnomad = f'https://gnomad.broadinstitute.org/region/__chrom__-__start__-__end__?dataset={gnomadGenome}'
    columns = list(df)
    linksUCSC = []
    linksGnomadRegion = []
    linksGnomadfwd = []
    linksGnomadrev = []

    for f,r,f_start,f_end,r_start,r_end,chrom in zip(df['fwd 5>3'],df['rev 5>3'], df['fwd start'].astype(str),df['fwd end'].astype(str),df['rev start'].astype(str),df['rev end'].astype(str),df['chrom'].astype(str)):
        link = l_ucsc
        link = link.replace('___frwdP___',f)
        link = link.replace('___reverseP___',r)
        linksUCSC.append(link)

        link = lgnomad
        link = link.replace('__chrom__',chrom)
        link = link.replace('__start__',f_start)
        link = link.replace('__end__',r_end)
        linksGnomadRegion.append(link)



        link = lgnomad
        link = link.replace('__chrom__',chrom)
        link = link.replace('__start__',f_start)
        link = link.replace('__end__',f_end)
        linksGnomadfwd.append(link)


        link = lgnomad
        link = link.replace('__chrom__',chrom)
        link = link.replace('__start__',r_start)
        link = link.replace('__end__',r_end)
        linksGnomadrev.append(link)

    df['UCSC PCR']= linksUCSC
    df['gnomAD region'] = linksGnomadRegion
    df['forward primer gnomad'] = linksGnomadfwd
    df['reverse primer gnomad'] = linksGnomadrev
    df['Nr. Abb'] = list(range(len(df)))
    fp = [n for n,x in enumerate(columns) if x =='fwd start'][0]
    rp = [n for n,x in enumerate(columns) if x =='rev start'][0]

    newC = ['Nr. Abb', 'chrom', 'transcript','gene', 'penalty', 'amplicon',  
    'fwd c. start', ' fwd c. end' ,'fwd 5>3' , 'fwd Tm', 'forward primer gnomad' ,
    ' rev c. start', ' rev c. end', 'rev 5>3', 'rev Tm',  'reverse primer gnomad',
    'UCSC PCR','gnomAD region', 
    'fwd len','fwd GC', 'fwd start','fwd end',
    'rev len','rev GC', 'rev start','rev end',
    'aln fwd 5>3', 'aln rev 5>3',
    'name','query']
    #newC = ['gnomAD region'] + ['UCSC PCR']+columns[:fp]+['forward primer gnomad'] + columns[fp:rp] + ['reverse primer gnomad'] + columns[rp:]
    #print(newC)
    df = df[newC]
    return df

