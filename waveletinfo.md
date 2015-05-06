### Approaches tried so far

#### Full discrete wavelet transform on the columns and rows of the image
Best results observed with db16 wavelets and periodic boundary conditions. Periodic boundary conditions cause issues when the image is non-sparse near the edges. Using an even extension of the image seems to be more appropriate, but implementation difficulties have impeded progress.

The locality gained from moving to a wavelet-based approach comes with the disadvantage of insensitivity to large scale influences. May be a fundamentally irreconcileable tradeoff.

#### Partial wavelet decomposition, with DCT used on the remaining approximate image
An attempt to get the large scale behavior of the DCT approach and the high frequency locality of the wavelet approach. Has not yet produced significant results.

#### 2D discrete wavelet transform into a lowpass approximation and horizontal, vertical, and diagonal detail images
This is the type of 2D wavelet transform implemented in Matlab's wavelet toolbox. Initial implementations using this type of transform showed very poor results.

#### DCT with a different norm
Changing the (3/2) exponent in the Lambda factor used to calculate the norm changes how much the norm attenuates high frequencies. An exponent of 1 increases locality to roughly match the best results from the wavelet approach, but sticking with the DCT avoids the implementation issues associated with wavelets.


