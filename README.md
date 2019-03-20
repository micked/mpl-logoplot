# mpl-logoplot
Simple utilities for logoplots in directly in Matplotlib

```python
fig, ax = plt.subplots()
mpl_logoplot.logoplot(ax, seq, 'ATGC', height=entropy, color='dna')
```