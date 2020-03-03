"""
Functions to count and cluster amino acid sequences.
"""

import numpy as np

from . import utils


class PSFM:
    """Meta class for a position specific scoring matrix"""
    def __init__(self, pssm, alphabet=utils.AMINO_ACIDS, comments=(), consensus=None):
        self._psfm = psfm
        self.alphabet = alphabet
        self.comments = comments
        self.consensus = consensus

    @classmethod
    def from_txt(cls, lines):
        """Load a scoring/frequency matrix from text"""
        alphabet = None
        dtype = int
        matlines = []
        comments = []
        consensus = []
        for line in lines:
            line = line.strip()
            if line.startswith('#'):
                comments.append(line.lstrip('#'))

            line = line.split('#')[0]
            if not line:
                continue

            if alphabet is None:
                alphabet = tuple(line.strip().split())
            else:
                if dtype is float or '.' in line:
                    dtype = float

                line = line.split()
                consensus.append(line[1])
                matlines.append([dtype(n) for n in line[2:]])

        mat = np.array(matlines)
        return cls(mat, alphabet=''.join(alphabet), comments=comments, consensus=''.join(consensus))

    def to_txt(self, mat=None):
        """Output a matrix in text format."""
        lines = []
        for line in self.comments:
            lines.append(f'#{line}')

        if mat is None:
            mat = self.pssm

        consensus = self.consensus
        if self.consensus is None:
            consensus = ['X'] * mat.shape[0]

        l = len(str(mat.shape[0]))
        lines.append(' ' * (l + 3) + ' '.join([f' {a:<5}' for a in self.alphabet]))
        for i in range(mat.shape[0]):
            line = f'{i+1:>{l}} {consensus[i]} ' + ' '.join([f'{n:>6.3f}' if n else ' 0    ' for n in mat[i]])
            lines.append(line)

        return lines

    def psfm(self):
        return self._psfm

    def pssm(self, bg='blosum62', return_psfm=False):
        if bg is None:
            bg = np.ones(len(self.alphabet)) / len(self.alphabet)
        elif bg in ('blosum62', 'blosum', 'bl62', 'bl'):
            bg = utils.bgfreq_array(self.alphabet)

        psfm = self.psfm()
        pssm = np.zeros_like(psfm, dtype='float')
        mask = psfm > 0
        pssm[mask] = np.log2((psfm / bg)[mask])

        if return_psfm:
            return pssm, psfm

        return pssm

    def shannon_logo(self):
        pssm, psfm = self.pssm(bg=None, return_psfm=True)
        return np.sum(pssm * psfm, axis=1, keepdims=True) * psfm

    def kullback_leibler_logo(self):
        pssm, psfm = self.pssm(bg='blosum62', return_psfm=True)
        return np.sum(pssm * psfm, axis=1, keepdims=True) * psfm * np.sign(pssm)

    def weighted_kullback_leibler_logo(self):
        pssm, psfm = self.pssm(bg='blosum62', return_psfm=True)
        return np.sum(pssm * psfm, axis=1, keepdims=True) * pssm / np.sum(np.abs(pssm), axis=1, keepdims=True)

    def p_weighted_kullback_leibler_logo(self):
        pssm, psfm = self.pssm(bg='blosum62', return_psfm=True)
        return np.sum(pssm * psfm, axis=1, keepdims=True) * pssm * psfm / np.sum(
            np.abs(pssm) * psfm, axis=1, keepdims=True)


class AlignmentPSFM(PSFM):
    """Position-specific scoring matrix from a set of sequences"""
    def __init__(self, sequences, alphabet=utils.AMINO_ACIDS, clustering='hobohm1', weight_on_prior=200, **kwargs):
        self.seqlen = None
        self.alphabet = str(alphabet)
        self.sequences = []

        #Validate sequence lengths
        for seq in sequences:
            seq = str(seq)
            self.sequences.append(seq)
            if self.seqlen == None:
                self.seqlen = len(seq)
            elif self.seqlen != len(seq):
                raise Exception('All sequences must be of same length!')

        if isinstance(clustering, str):
            if clustering.lower() in ('hobohm', 'hobohm1'):
                clustering = hobohm1_factory(kwargs.get('hobohm1_threshold', 0.63))
            elif clustering.lower() in ('heuristic', ):
                clustering = heuristic
            else:
                raise Exception(f'Unknown clustering method: {clustering}')

        self.clustering = clustering
        self.weight_on_prior = weight_on_prior
        self.consensus = None
        self.comments = []
        self.gaps = '-.'

    @property
    def alphabet_index(self):
        alpha_index = {a: i for i, a in enumerate(self.alphabet, 1)}

        for gap in self.gaps:
            alpha_index[gap] = 0

        return alpha_index

    def alignment_array(self):
        if set(self.gaps).intersection(set(self.alphabet)):
            raise Exception(f'Alphabet ({self.alphabet}) must not contain gaps ({self.gaps})')

        alphabet_index = self.alphabet_index
        alignment_arr = []
        for seq in self.sequences:
            seq_coded = []
            for letter in seq:
                seq_coded.append(alphabet_index[letter])
            alignment_arr.append(seq_coded)

        return np.array(alignment_arr)

    def sequence_weights(self):
        if self.clustering is None:
            return np.ones(len(self.sequences))

        return self.clustering(self.alignment_array())

    def psfm(self):
        """Position specific frequency matrix"""
        a_index = self.alphabet_index
        sequence_weights = self.sequence_weights()

        countmat = np.zeros([len(self.sequences[0]), len(self.alphabet)], dtype='float')

        for j, seq in enumerate(self.sequences):
            for i, a in enumerate(seq):
                a_idx = a_index[a] - 1
                if a_idx >= 0:
                    countmat[i, a_idx] += sequence_weights[j]

        frequency_matrix = countmat / np.sum(sequence_weights)

        if self.weight_on_prior:
            blosum_mat = utils.blosum62_array(self.alphabet)
            pseudo_counts = np.dot(frequency_matrix, blosum_mat)
            a = np.sum(sequence_weights) - 1
            b = self.weight_on_prior
            frequency_matrix = (a * frequency_matrix + b * pseudo_counts) / (a + b)

        return frequency_matrix


#
# Clustering/weighting


def hobohm1_factory(threshold):
    """Return Hobohm1 clustering func with a given threshold"""
    def _hobohm1(alignment_array):
        return hobohm1(alignment_array, threshold)

    return _hobohm1


def hobohm1(alignment_array, threshold=0.63):
    """Hobohm1 as implemented in seq2logo"""
    n_seqs, seqlen = alignment_array.shape
    seqsort = np.argsort(np.sum(alignment_array == 0, axis=1))
    alignment_array = alignment_array[seqsort]

    clusters = np.arange(n_seqs)
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            if clusters[j] != j:
                continue
            sim = alignment_array[i] == alignment_array[j]
            sim[alignment_array[i] == 0] = False
            sim = np.sum(sim) / np.sum(alignment_array[i] != 0)
            if sim >= threshold:
                clusters[j] = i

    counts = np.bincount(clusters)
    weights = np.ones(n_seqs) / counts[clusters]

    seqsort_undo = np.empty(seqsort.size, 'int')
    seqsort_undo[seqsort] = np.arange(seqsort.size)
    return weights[seqsort_undo]


def heuristic(alignment_array):
    """Heuristic (position-based) clustering.
    
    Reference: https://doi.org/10.1016/0022-2836(94)90032-9
    
    """
    n_seqs, seqlen = alignment_array.shape

    #Calculate weight per position
    #weight(seq, pos) = 1 / (count(letter_seq)_pos * n_uniq_letters_pos)
    weight = np.zeros_like(alignment_array, dtype='float')
    for pos in range(seqlen):
        uniq, indices, counts = np.unique(alignment_array[:, pos], return_counts=True, return_inverse=True)
        weight[:, pos] = (1 / (counts * uniq.size))[indices]

    return np.sum(weight, axis=1)
