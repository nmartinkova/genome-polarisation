# Genome polarisation with diemr

This repository hosts the source files for the Bookdown guide  
**“Genome polarisation with diemr”**  
published at https://nmartinkova.github.io/genome-polarisation/.

---

## About the book

The book provides a practical guide to genome polarisation using the  
[`diemr`](https://github.com/nmartinkova/diemr) R package.  
It explains how to:
- import and polarise genotypes,
- filter markers by diagnostic index,
- smooth polarised genomic states along chromosomes,
- compute hybrid indices and visualise genomic structure.

Each chapter includes example code and interpretation notes based on the  
*diagnostic index expectation maximisation (diem)* algorithm  
described in  
[Baird *et al.* (2023)](https://doi.org/10.1111/2041-210X.14010).

---

## Repository contents

- `.Rmd` files — source text for each chapter  
- `_bookdown.yml`, `_output.yml` — build configuration  
- `docs/` — rendered HTML files for GitHub Pages  
- `figs/` — pre-rendered figures included in the book  

To render locally:

```r
bookdown::render_book("index.Rmd")
```