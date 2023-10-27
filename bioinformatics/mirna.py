import re
import pandas as pd
from Bio import SeqIO
import numpy as np

class MicroRNA:
    """
    MicroRNA object holding structure information for named miRNA
    
    Args
    structure    structure formatted as txt
    name         unique miRNA name
    """
    
    def __init__(self, structure, name):
        self.structure = structure
        self.name = name
        
    def create_fasta_matures(self):
        fivep, threep, _ = self.extract_mature()
        return f'>{self.name}\n{fivep}&{threep}'
    
    def create_fasta_hairpin(self):
        fivep, threep, _ = self.extract_hairpin()
        return f'>{self.name}\n{fivep}&{threep}'
        
    def extract_true_bpp(self, df, seqs = None, bp = None, mature = False, extract = 'max'):
        """
        Function to create x, y array of basepair binding probabilities. 
        x < 0: 5' strand, x > 0: 3' strand
        Defaults to maximum binding probability if no bond in structure.

        Args:
        df      Pandas df holding basepair binding probabilities
        seqs    Dictionary containing sequence of 5', 3' 
        bp      True basepair bonds from structure. Defults None
        mature  Bool, extract mature sequences from structure. Default False
        extract String, ['max', 'min', 'mean', 'median', 'sum'], method of deciding 
                probability in absence of true bond

        Returns:
        x       numpy array of x coordinates for plotting
        y       numpy array of basepair binding probability
        """

        #extract unique indices for each strand
        unique_5p = df['i'].unique()
        unique_3p = df['j'].unique()
        
        if seqs == None:
            if mature:
                fivep, threep, bp = self.extract_mature()
            else:
                fivep, threep, bp = self.extract_hairpin()
            seqs = {'5p': fivep, '3p': threep}

        #make empty array, x
        y = np.zeros(len(seqs['5p']) + len(seqs['3p']) + 1)
        x = np.linspace(-(len(seqs['5p'])), len(seqs['3p']), len(y))

        #fill y with basepair binding probabilities
        for i in range(len(unique_5p)):
            yx = int(unique_5p[i])-1
            a = df.loc[df['i'] == unique_5p[i]]
            tups = [(int(a['i'].iloc[k]), int(a['j'].iloc[k])) for k in range(len(a['i']))]
            for j in range(len(tups)):
                if tups[j] in list(bp):
                    y[yx] = float(a['prob'].iloc[j])
                    break
                else:
                    b = pd.Series(a['prob'], dtype = 'float')
                    y[yx] = self._extracts(b, extract)
                     

        for i in range(len(unique_3p)):
            yx = int(unique_3p[i])
            a = df.loc[df['j'] == unique_3p[i]]
            tups = [(int(a['i'].iloc[k]), int(a['j'].iloc[k])) for k in range(len(a['i']))]
            for j in range(len(tups)):
                if tups[j] in list(bp):
                    y[yx] = float(a['prob'].iloc[j])
                    break
                else:
                    b = pd.Series(a['prob'], dtype = 'float')
                    y[yx] = self._extracts(b, extract)
                     
                    
        return x, y
    
    def extract_max_bpp(self, df, seqs = None, mature = False, extract = 'max'):
        """
        Function to create x, y array of basepair binding probabilities. 
        x < 0: 5' strand, x > 0: 3' strand
        Maximum binding probability with no regard to same base binding to multiple 
        other bases

        Args:
        df      Pandas df holding basepair binding probabilities
        seqs    Dictionary containing sequence of 5', 3' 
        mature  Bool, extract mature sequences from structure. Default False
        extract String, ['max', 'min', 'mean', 'median', 'sum'], method of deciding 
                probability in absence of true bond

        Returns:
        x       numpy array of x coordinates for plotting
        y       numpy array of basepair binding probability
        """

        #extract unique indices for each strand
        unique_5p = df['i'].unique()
        unique_3p = df['j'].unique()
        
        if seqs == None:
            if mature:
                fivep, threep, _ = self.extract_mature()
            else:
                fivep, threep, _ = self.extract_hairpin()
            seqs = {'5p': fivep, '3p': threep}

        #make empty array, x
        y = np.zeros(len(seqs['5p']) + len(seqs['3p']) + 1)
        x = np.linspace(-(len(seqs['5p'])), len(seqs['3p']), len(y))
        
        #fill y with basepair binding probabilities
        for i in range(len(unique_5p)):
            yx = int(unique_5p[i]) - 1
            a = pd.Series(df['prob'].loc[df['i'] == unique_5p[i]], dtype = 'float')
            y[yx] = self._extracts(a, extract)

        for i in range(len(unique_3p)):
            yx = int(unique_3p[i])
            a = pd.Series(df['prob'].loc[df['j'] == unique_3p[i]], dtype = 'float')
            y[yx] = self._extracts(a, extract)

        return x, y
    
    def _extracts(self, series, extract):
        """
        function to extract bpp
        """
        if extract == 'max':
            return series.max()
        elif extract == 'min':
            return series.min()
        elif extract == 'mean':
            return series.mean()
        elif extract == 'median':
            return series.median()
        elif extract == 'sum':
            return np.sum(np.array(series))
    
    
    def extract_mature(self):
        """
        Function to extract mature sequence from string structure pre-miRNA,
        record match/non-match
        """

        split = re.split(r'\n', self.structure)

        mm = []

        fivep = ''
        threep = ''

        #find index of last section of stem region
        i = -1
        while i >= -len(split[2]):
            if split[2][i] == '|':
                ix = i
                break
            i -= 1

        #setup 3' indexing
        j = 2*(len(split[0])+ix+1)
        #5' strand
        for i in range(len(split[0])+ix+1):
            if split[0][i] == ' ':
                if split[1][i].isupper():
                    fivep += split[1][i]
                    mm.append([i,j])
                else:
                    continue
            elif split[0][i] == '-':
                continue
            elif split[0][i].isupper():
                fivep += split[0][i]
                mm.append(-1)
        #3' strand
        i = len(split[0]) + ix
        while i >= 0:
            if split[-1][i] == ' ':
                if split[-2][i].isupper():
                    threep += split[-2][i]
                else:
                    i -= 1
                    continue
            elif split[-1][i] == '-':
                i -= 1
                continue
            elif split[-1][0].isupper():
                threep += split[-1][i]
            i -= 1


        #re-index mm
        for i in range(len(mm)):
            if mm[i] == -1:
                continue
            else:
                if len(fivep) + len(threep) - i -1 != len(fivep):
                    mm[i][0] = i
                    mm[i][1] = len(fivep) + len(threep) - i -1
                    mm[i] = tuple(mm[i])
                else:
                    mm[i] = -1
        for i in range(len(mm)):
            if mm[i] != -1:
                if mm[i][0] == mm[i][1]:
                    mm[i] = -1

        return fivep, threep, mm
    
    def extract_hairpin(self):
        """
        Function to extract mature sequence from string structure pre-miRNA,
        record match/non-match
        """

        split = re.split(r'\n', self.structure)

        mm = []

        fivep = ''
        threep = ''

        #find index of last section of stem region
        i = -1
        while i >= -len(split[2]):
            if split[2][i] == '|':
                ix = i
                break
            i -= 1

        #setup 3' indexing
        j = 2*(len(split[0])+ix+1)
        #5' strand
        for i in range(len(split[0])+ix+1):
            if split[0][i] == ' ':
                fivep += split[1][i]
                if split[2][i] == '|':
                    mm.append([i,j])
                else:
                    mm.append(-1)
            elif split[0][i] == '-':
                continue
            else:
                fivep += split[0][i]
                if split[2][i] == '|':
                    mm.append([i,j])
                else:
                    mm.append(-1)
        #3' strand
        i = len(split[0]) + ix
        while i >= 0:
            if split[-1][i] == ' ':
                threep += split[-2][i]
            elif split[-1][i] == '-':
                i -= 1
                continue
            else:
                threep += split[-1][i]
            i -= 1


        #re-index mm
        for i in range(len(mm)):
            if mm[i] == -1:
                continue
            else:
                if len(fivep) + len(threep) - i -1 != len(fivep):
                    mm[i][0] = i
                    mm[i][1] = len(fivep) + len(threep) - i -1
                    mm[i] = tuple(mm[i])
                else:
                    mm[i] = -1
        for i in range(len(mm)):
            if mm[i] != -1:
                if mm[i][0] == mm[i][1]:
                    mm[i] = -1

        return fivep, threep, mm
    
    
def parse_postscript(ps, name):
    """
    Function to parse post-script file to extract base pair binding probabilities from RNAcofold
    
    Args:
    ps:   post-script file
    name: name to save .csv to, str (format: <path>/<name>.csv
    """
    
    #use regex to isolate the square roots of probabilities
    t2 = re.search(r'%start of base pair probability data\n(.*\n*)*\nshowpage\nend', ps, re.S)
    t2 = re.sub(r'%!PS.*?start of base pair probability data\n2', '', t2.group(), flags = re.S)
    t2 = re.split(r'\n', t2)
    t2 = t2[2:-3]
    
    expt = [re.split(r' ', i) for i in t2]
    
    out = {'i': [],
        'j': [],
        'prob': []}
    
    for row in expt:
        out['i'].append(row[0])
        out['j'].append(row[1])
        out['prob'].append(row[2])
    
    df = pd.DataFrame(out)
    df.to_csv(name, index = False)
    
def calc_stats(data_dict, subset = None):
    """
    Function to calculate means, sem for all mirnas in data_dict
    
    data_dict format:
    data = {mir{'x': [...], 'y': [...]}}
    
    Args:
    data_dict  Dictionary of x, y for all mirs
    subset     List of mirnames to subset from data_dict
    
    Returns:
    mean    list, bpp
    sem     list, bpp
    x       list, bp position
    
    """
    import numpy as np
    import scipy.stats as st
    
    mirs = _extract_data(data_dict, subset = subset)

    x = np.array(list(mirs.keys()), dtype = 'int')
    x.sort()
    
    y_mean = []
    y_sem = []
    for ix in x:
        y_mean.append(np.mean(np.array(mirs[str(ix)])))
        y_sem.append(st.sem(np.array(mirs[str(ix)])))
        
    return y_mean, y_sem, x

    
def calc_auc(data_dict, subset = None):
    """
    Calculate area under curve for miRNAs using Simpson's rule
    Normalises for miRNA length
    """
    import numpy as np
    from scipy.integrate import simpson
    
    auc = {}
    for mir in data_dict:
        _dat = data_dict[mir]
        
        if subset != None:
            if mir.lower() in [a.lower() for a in subset]:
                x = np.array(list(_dat['x']), dtype = 'int')
                x.sort()
                y = []
                df_temp = pd.DataFrame(_dat)
                for ix in x:
                    y.append(float(df_temp['y'].loc[df_temp['x'] == ix]))
                
                auc[mir] = simpson(y = y,
                               x = x,
                               even = 'avg')/len(x)
            
        else:
            x = np.array(list(_dat['x']), dtype = 'int')
            x.sort()
            y = []
            df_temp = pd.DataFrame(_dat)
            for ix in x:
                y.append(float(df_temp['y'].loc[df_temp['x'] == ix]))
        
            auc[mir] = simpson(y = y,
                               x = x,
                               even = 'avg')/len(x)
        
    return auc

def _extract_data(data_dict, subset = None):
    mirs = {}
    
    for mir in data_dict:
        _dat = data_dict[mir]
        for i in range(len(_dat['x'])):
            keys = mirs.keys()
                           
            if subset != None:
                subset = list(subset)
                if mir.lower() in [a.lower() for a in subset]:
                    
                    if str(int(_dat['x'][i])) not in keys:
                           mirs[str(int(_dat['x'][i]))] = []
                    mirs[str(int(_dat['x'][i]))].append(_dat['y'][i])
            
            else:
                if str(int(_dat['x'][i])) not in keys:
                    mirs[str(int(_dat['x'][i]))] = []
                mirs[str(int(_dat['x'][i]))].append(_dat['y'][i])
    
    return mirs