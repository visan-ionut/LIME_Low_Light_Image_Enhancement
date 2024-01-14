================================
LIME-Low Light Image Enhancement
================================

The proposed method employs a straightforward yet effective
approach by individually estimating the illumination of each
pixel. The illumination map is initially computed by finding
the maximum intensity in the R, G, and B channels, followed
by refinement with a structured approach. This refined illumination
map serves as a foundation for enhancing low-light images,
showcasing superior results compared to several state-of-the-art methods.
The project falls under the Retinex-based category, utilizing Multiple
Augmented Lagrange (ALM) algorithms to precisely refine the illumination
map and reduce computational costs. The implementation is versatile,
demonstrating adaptability to various structural weighting strategies,
making it a valuable tool for enhancing image quality in low-light scenarios.
