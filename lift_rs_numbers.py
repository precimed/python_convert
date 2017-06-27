import gzip

class LiftRsNumbers:
    def _read_rs_history(self, histFile):
        RS_HISTORY = set() # store rs

        print("Reading '{}' file...".format(histFile))
        for ln in gzip.open(histFile, mode='rt'):
            fd = ln.strip().split('\t')
            # Some very few entries in SNPHistory file are about
            # re-activation SNPs (not about deleting them).
            # We just need to ignore those entries
            if ln.lower().find('re-activ') < 0:
                RS_HISTORY.add(fd[0])
        print('{} entries found'.format(len(RS_HISTORY)))
        return RS_HISTORY

    def _read_rs_merge(self, mergFile):
        RS_MERGE = dict() # high_rs -> (lower_rs, current_rs)

        print("Reading '{}' file...".format(mergFile))
        for ln in gzip.open(mergFile, mode='rt'):
            fd = ln.strip().split('\t')
            h, l = fd[0], fd[1]
            c = fd[6]
            RS_MERGE[h] = (l, c)

        print('{} entries found'.format(len(RS_MERGE)))
        return RS_MERGE

    def __init__(self, hist_file=None, merge_file=None):
        self._RS_HISTORY = self._read_rs_history(hist_file)
        self._RS_MERGE = self._read_rs_merge(merge_file)

    def lift(self, rsvec):
        unchanged = 0; lifted = 0; deleted = 0; not_rs_number = 0;
        RS_LIFTED = rsvec.copy(); nsnps = len(rsvec)
        print("Lifting rs# numbers for n={} SNPs...".format(nsnps))
        is_rs_number = [x.startswith('rs') and x[2:].isdigit() for x in rsvec]
        rsvec = [x[2:] for x in rsvec]
        for i in range(nsnps):
            if not is_rs_number[i]:
                not_rs_number += 1
                continue
            rs = rsvec[i]
            if rs not in self._RS_MERGE and rs not in self._RS_HISTORY:
                unchanged += 1
                continue
            while True:
                if rs in self._RS_MERGE:
                    rsLow, rsCurrent = self._RS_MERGE[rs]
                    if rsCurrent not in self._RS_HISTORY and rsCurrent != '':
                        RS_LIFTED[i] = 'rs' + rsCurrent; lifted += 1
                        break
                    else:
                        rs = rsLow
                else:
                    # Such SNPs were deleted from SNPdb,
                    # look it up here: https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=71800898
                    RS_LIFTED[i] = None; deleted += 1
                    break

        return RS_LIFTED, {'invalid rs#':not_rs_number, 'unchanged':unchanged, 'lifted':lifted, 'deleted':deleted}
