The import is fixed, and the notebook already covers a lot: sample/count validation, library sizes, zero counts, low-expression filtering, VST, and sample correlation checks. Based on the notebook currently in the repository, I’d address these gaps before calling QC complete:

- Add PCA on the VST matrix. Color samples by Condition and, if available, batch, tissue type, or sequencing metadata. This can show whether samples cluster by condition or by a technical factor.

- Review sample outliers. Use the PCA and correlation plots to identify unusual samples. Record any exclusions and the reason; don’t remove a sample just because it looks different.

- Check batch effects if you have batch metadata. The current metadata has only Condition, so the notebook cannot assess whether a technical batch is driving patterns. If that information is unavailable, note that limitation.

- Verify what “Healthy” means. Document whether those samples are independent healthy donors or adjacent normal tissue, since that changes how you should describe the comparison.

- Clarify the scope of filtering. The current rule—at least 10 counts in at least 10% of samples—is a reasonable exploratory filter. For future cross-validated ML, though, do feature selection and other data-learned preprocessing inside each training fold to avoid leakage.

- Clean up unused code and notebook state. save_figure is imported but not used, and several outputs appear blank in the saved notebook. After settling the analysis, rerun it from top to bottom and save the outputs so the notebook shows the results you reviewed.


I’d consider PCA, an outlier review, and a written note about batch information and sample type the main remaining QC work. Differential expression analysis is a separate next analysis, useful for finding biologically associated genes, but it isn’t required to finish this QC notebook or start a carefully validated classification baseline.