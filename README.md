# MXRGA

This project demonstrates how to utilize the PIPER docking program along with symmetry handling implemented in MATLAB to predict the protein crystal packing of a given protein. It also exhibits how to assess potential conditions for protein crystallization by analyzing the similarity of protein packing interfaces calculated using AAI-PatchBag.

## Contents:

- ### Case study
  #### * Please follow the steps according to the folder numbers, and pay attention to the file directories.
  (1) P41212 and P43212 samples using MXRGA (from sample download to simulated crystal packing models)
     - Input:  PDB file or AF2 protein structure of a given sample <br>
     - Output: Simulated crystal packing models                    <br>
  
  (2) Protein crystallization condition ranking through the similarity of protein packing interfaces evaluated by AAI-PatchBag
     - Input:  PDB_ID (Please rename it in Script.m) and corresponding PDB file <br>
     - Output: AAI-PatchBag based condition (PDB ID) ranking in mat file        <br>

- ### Methodology
  (1) Single-chain human protein download and sample reduction <br>
  (2) Protein crystallization condition extraction             <br>
  (3) Crystallization condition preprocessing                  <br>
  (4) DIPER coefficients and PackQ training                    <br>
  (5) Establish of AAI-PatchBag                                <br>
  (6) Correlation of PDB crystallization condition and protein packing  <br>

## Requirements:

1. PIPER    : physical-based docking program.                                <br>
   available in: https://cluspro.org/downloads.php
2. EDTSurf  : rapid macromolecular surface calculation.                      <br>
   available in: https://zhanggroup.org/EDTSurf/
3. generate_3dzd.py  : molecular shape in ply format and 3DZD computation.   <br>
   available in: https://github.com/kiharalab/3d-af_surfer
4. MakeShape, Shape2Zernike	: 3DZD computation.                              <br>
   available in: https://github.com/jerhoud/zernike3d/tree/main
