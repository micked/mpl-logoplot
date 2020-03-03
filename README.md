# mpl-logoplot
Simple script for logoplots directly in Matplotlib

```python
fig, ax = plt.subplots()
mpl_logoplot.logoplot(ax, seq, 'ATGC', height=entropy, color='dna')
```

![Random DNA](https://raw.githubusercontent.com/micked/mpl-logoplot/master/Examples/DNA.png)


Experimental support for generating position specific frequency/scoring matrices (PSFM/PSSM):

```python
alignment = mpl_logoplot.psfm.AlignmentPSFM(seqs, alphabet=mpl_logoplot.AMINO_ACIDS, clustering='hobohm', weight_on_prior=200)

logo = alignment.shannon_logo()
mpl_logoplot.logoplot(ax, logo, alignment.alphabet, color='protein')
```

![Peptide](https://raw.githubusercontent.com/micked/mpl-logoplot/master/Examples/peptide.png)