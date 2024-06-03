# SSCTV model Performance

## HSI Denoising Benchmark

Our SSCTVmodel has been rigorously tested in the task of sparsity noise reduction on HSI (Hyperspectral Imaging) examples. The empirical results demonstrate a significant improvement over the baseline models. We have benchmarked the performance using widely accepted metrics such as PSNR (Peak Signal-to-Noise Ratio), SSIM (Structural Similarity Index), and ERGAS (Erreur Relative Globale Adimensionnelle de Synthèse).

Run `Simulation` for the experiments with the simulated data we generated, and for denoising experiments on real data, please run two demos.

After running `Demo_of_HSI_Denoising`, our model achieved the following metrics:

- **PSNR**: Our model exhibits superior noise reduction capabilities, resulting in higher mPSNR values compared to the baseline models.
- **SSIM**: The structural integrity of the denoised images is maintained more effectively by our model, as reflected by its higher mSSIM scores.
- **ERGAS**: The lower ERGAS values obtained from our model indicate a smaller relative global error in synthesis, implying a more accurate reconstruction of the original image.

What's more, we supply more experiment about multi spectral image in CAVE dataset by running `Demo_of_MSI_Denoising`，please download it on https://www1.cs.columbia.edu/CAVE/databases/multispectral/.

## Anomaly Detection Performance

In the domain of Anomaly Detection, we have also assessed our model's performance by computing the Area Under the ROC Curve (AUC). Our model's AUC scores are indicative of its robustness and precision in identifying anomalies when compared to baseline models. The higher AUC value underscores the enhanced detection capabilities of our SSCTVmodel.

The comparative results highlight the efficacy of the SSCTVmodel in providing state-of-the-art performance for both HSI denoising and anomaly detection tasks.

We welcome you to test our model with your datasets and compare its performance against other models.



