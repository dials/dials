``dials.find_spots``: added ``radial_profile`` threshold algorithm. This
calculates an average background in 2θ shells and identifies peak pixels
at a user-controllable level above the background. This simple method
is particularly appropriate for cases with strong rotationally-symmetric
background, such as electron diffraction images. An optional blurring
function helps to suppress noise peaks and to join split spots.
