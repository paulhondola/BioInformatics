# Obiectiv

Această temă extinde **clustering-ul** cu analiza **filogenetică** realizată în Lab 4.  
Veți re-folosi dataset-ul propriu și arborele filogenetic creat anterior, pentru a compara clusterele statistice cu ramurile evolutive.  

---

## Instrucțiuni

### Pasul 1 — Dataset și Arbore Filogenetic
- Refolosiți **dataset-ul multi-FASTA** construit în Lab 4.  
- Folosiți arborele filogenetic salvat (`tree_<handle>.nwk`) ca referință.  

### Pasul 2 — Clustering
- Preprocesați datele: transformați-le într-o **matrice de distanțe / similaritate** sau vectori de caracteristici.  
- Aplicați cel puțin **două metode** de clustering:  
  - Hierarchical (average linkage sau altă metodă la alegere)  
  - K-means (testați mai multe valori pentru K)  
  - DBSCAN (variați eps / min_samples)  
- Evaluați calitatea clusterelor folosind un indice (ex: Silhouette Score).  

### Pasul 3 — Integrare cu Filogenetica
- Comparați clusterele cu ramurile arborelui filogenetic.  
- Analizați dacă genele/ secvențele din același cluster apar și în același clade.  
- Discutați diferențele observate:  
  - suprapuneri,  
  - discrepanțe,  
  - ipoteze biologice (ex: evoluție convergentă, diversificare funcțională).  

---

## Livrabile
1. **Codul Python** folosit pentru clustering și comparație cu arborele.  
2. **Vizualizări**:  
   - dendrogramă,  
   - ploturi PCA pentru clustering,  
   - arbore filogenetic adnotat cu clustere.  
3. **Raport scurt (max 2 pagini)** care să includă:  
   - Rezultatele clustering-ului,  
   - Compararea cu arborele filogenetic,  
   - Interpretarea biologică/funcțională (vezi intrebari de reflectie).  

---

## Notare
- **Dataset și arbore reutilizat corect ** — 1p  
- **Aplicarea corectă a două metode de clustering** — 3p  
- **Vizualizări clare (dendrogramă, PCA, arbore adnotat)** — 3p  
- **Analiza comparativă clustering vs. filogenetică** — 2p  
- **Calitatea raportului (structură, claritate, reflecție critică)** — 1p  
- **Bonus (+1p):** Folosiți 3 metode de clustering și discutați diferențele între ele.  

---

## Întrebări de reflecție
- Cum se aliniază rezultatele clustering-ului cu arborele filogenetic?  
- Ce explică diferențele observate între metode și între tree vs. clustering?  
- Cum poate fi utilă combinarea clustering-ului cu filogenetica în alte aplicații (ex: identificarea subtipurilor de boală, gene families, evoluție funcțională)?
