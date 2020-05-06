% minFunc
fprintf('Compiling minFunc files...\n');
mex -compatibleArrayDims -outdir minFunc minFunc/mcholC.c
mex -compatibleArrayDims -outdir minFunc minFunc/lbfgsC.c
mex -compatibleArrayDims -outdir minFunc minFunc/lbfgsAddC.c
mex -compatibleArrayDims -outdir minFunc minFunc/lbfgsProdC.c

% KPM
fprintf('Compiling KPM files...\n');
mex -compatibleArrayDims -IKPM -outdir KPM KPM/max_mult.c

% DAGlearn
fprintf('Compiling DAGlearn files...\n');
mex -compatibleArrayDims -IDAGlearn/ancestorMatrix -outdir DAGlearn/ancestorMatrix DAGlearn/ancestorMatrix/ancestorMatrixAddC_InPlace.c
mex -compatibleArrayDims -IDAGlearn/ancestorMatrix -outdir DAGlearn/ancestorMatrix DAGlearn/ancestorMatrix/ancestorMatrixBuildC.c

% L1GeneralOverlapping Group
fprintf('Compiling L1GeneralGroup files...\n');
mex -compatibleArrayDims -outdir L1GeneralGroup/mex L1GeneralGroup/mex/projectRandom2C.c
mex -compatibleArrayDims -outdir L1GeneralGroup/mex L1GeneralGroup/mex/auxGroupLinfProjectC.c
mex -compatibleArrayDims -outdir L1GeneralGroup/mex L1GeneralGroup/mex/auxGroupL2ProjectC.c

% L1GeneralOverlapping Group
fprintf('Compiling L1GeneralOverlapping Group files...\n');
mex -compatibleArrayDims -outdir L1GeneralOverlappingGroup L1GeneralOverlappingGroup/projectNDgroup_DykstraFastC.c

% UGM
fprintf('Compiling UGM files...\n');
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeEdgeVEC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ExactC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ExactC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_ChainC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeClampedPotentialsC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ICMC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_GraphCutC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Sample_GibbsC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_MFC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_LBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_LBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Infer_TRBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_TRBPC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_makePotentialsC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_PseudoNLLC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_LogConfigurationPotentialC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_AlphaExpansionC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_AlphaExpansionBetaShrinkC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_CRF_NLLC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_Decode_ChainC.c
mex -compatibleArrayDims -IUGM/mex -outdir UGM/compiled UGM/mex/UGM_makeCRFmapsC.c

% LLM2
fprintf('Compiling LLM2 files...\n');
mex -compatibleArrayDims -outdir LLM2/mex LLM2/mex/LLM2_inferC.c
mex -compatibleArrayDims -outdir LLM2/mex LLM2/mex/LLM2_suffStatC.c
mex -compatibleArrayDims -outdir LLM2/mex LLM2/mex/LLM2_pseudoC.c

% LLM
fprintf('Compiling LLM files...\n');
mex -compatibleArrayDims -outdir LLM/mex LLM/mex/LLM_inferC.c
mex -compatibleArrayDims -outdir LLM/mex LLM/mex/LLM_suffStatC.c
mex -compatibleArrayDims -outdir LLM/mex LLM/mex/LLM_pseudoC.c

% UGMep
fprintf('Compiling UGMep files\n');
mex -compatibleArrayDims -outdir UGMep/mex UGMep/mex/UGMep_Decode_ICMC.c
mex -compatibleArrayDims -outdir UGMep/mex UGMep/mex/UGMep_EnergyC.c
mex -compatibleArrayDims -outdir UGMep/mex UGMep/mex/UGMep_makeClampedEnergyC.c
mex -compatibleArrayDims -outdir UGMep/mex UGMep/mex/UGMep_Decode_GraphCutC.c
mex -compatibleArrayDims -outdir UGMep/mex UGMep/mex/UGMep_Decode_AlphaExpansionC.c
mex -compatibleArrayDims -outdir UGMep/mex UGMep/mex/UGMep_Decode_ExpandShrinkC.c

% Ewout
fprintf('Compiling Ewout files...\n');
mex -compatibleArrayDims -IForeign/Ewout -outdir Foreign/Ewout Foreign/Ewout/projectBlockL1.c Foreign/Ewout/oneProjectorCore.c Foreign/Ewout/heap.c
mex -compatibleArrayDims -IForeign/Ewout -outdir Foreign/Ewout Foreign/Ewout/projectBlockL2.c

% Misc
mex -compatibleArrayDims -outdir misc misc/sampleDiscrete_cumsumC.c

% SAG
mex -compatibleArrayDims -outdir SAG/mex SAG/mex/SGD_logistic.c -largeArrayDims
mex -compatibleArrayDims -outdir SAG/mex SAG/mex/ASGD_logistic.c -largeArrayDims
mex -compatibleArrayDims -outdir SAG/mex SAG/mex/PCD_logistic.c -largeArrayDims
mex -compatibleArrayDims -outdir SAG/mex SAG/mex/DCA_logistic.c -largeArrayDims
mex -compatibleArrayDims -outdir SAG/mex SAG/mex/SAG_logistic.c -largeArrayDims
mex -compatibleArrayDims -outdir SAG/mex SAG/mex/SAGlineSearch_logistic.c -largeArrayDims
mex -compatibleArrayDims -outdir SAG/mex SAG/mex/SAG_LipschitzLS_logistic.c -largeArrayDims
%mex -compatibleArrayDims -outdir SAG/mex SAG/mex/SGD_logistic_BLAS.c -largeArrayDims -lmwblas
%mex -compatibleArrayDims -outdir SAG/mex SAG/mex/ASGD_logistic_BLAS.c -largeArrayDims -lmwblas
%mex -compatibleArrayDims -outdir SAG/mex SAG/mex/PCD_logistic_BLAS.c -largeArrayDims -lmwblas
%mex -compatibleArrayDims -outdir SAG/mex SAG/mex/DCA_logistic_BLAS.c -largeArrayDims -lmwblas
%mex -compatibleArrayDims -outdir SAG/mex SAG/mex/SAG_logistic_BLAS.c -largeArrayDims -lmwblas
%mex -compatibleArrayDims -outdir SAG/mex SAG/mex/SAGlineSearch_logistic_BLAS.c -largeArrayDims -lmwblas
%mex -compatibleArrayDims outdir SAG/mex SAG/mex/SAG_LipschitzLS_logistic_BLAS.c -largeArrayDims -lmwblas

