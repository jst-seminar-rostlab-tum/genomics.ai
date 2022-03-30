# jst-genomics.ai
An open source web application to annotate and visualise single-cell sequencing data using neural networks (e.g. for cancer research) 



# General info
Genomics.ai is a platform to annotate and visualise single-cell sequencing data using neural networks (e.g. for cancer research).
The visualisation of annotation results is realised using a UMAP. 


# Visualization Process
1. Upload single-cell annotation data for visualisation
2. Process data using ML model on backend
3. Send csv file for visualisation to use
4. visualisation of annotation is generated using D3.js


# Main features currently implemented
- Create profile using email, (currently restricted to academic emails)
- Multipart file upload in chunks (necessary for larger files)
- Create multiple visualisation projects
- Processing of input using ML model
- Access to previous projects
- Shareable link of visualisation
- Classification of visualisation job ranging from “uploading” to “complete”.
