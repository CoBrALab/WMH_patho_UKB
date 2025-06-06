# Description

Analysis code for the manuscript "Characterizing spatiotemporal white matter hyperintensity pathophysiology in vivo to disentangle vascular and neurodegenerative contributions"

We share all analysis code, code used for generating plots, group-level raw results (e.g., 3D volume maps for spatial clusters, linear model outputs, etc.) and visualizations included in the paper.

No individual-level data from UKB or ADNI can be shared for confidentiality reasons.

For questions/comments, please reach out to Olivier Parent (olivier.parent@mail.mcgill.ca)

# Data available

Description of the files made available in data/. Files are available for download as a GitHub Release. All files are shared in both minc and nifti format.

Custom UK Biobank template space (T1w)
- UKB_template: 1mm isotropic
- UKB_template_2mm: 2mm isotropic
- UKB_template_mask: mask

Final data-driven WMH parcellation based on pathophysiology
- WMH_parc_patho: in UKB space
- WMH_parc_patho_MNI: in MNI ICBM152 non-linear symmetric 09c space

Prevalence maps for WMH and NAWM
- NAWM_prevalence: normal-appearing white matter labels prevalence
- WMH_prevalence: white matter hyperintensity labels prevalence
- WMH_mask: mask of voxels with >1 WMH label

Other WMH parcellation used for comparisons (all in UKB space)
- WMH_parc_vascular: parcellation of arterial territories (from Liu et al., 2023, Scientific Data, https://doi.10.1038/s41597-022-01923-0)
- WMH_parc_lobar: lobar parcellation
- WMH_parc_pv_deep: periventricular/deep parcellation

