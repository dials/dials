from __future__ import annotations

import random

from dials.array_family import flex


class Simulator:
    """Class to help with simulation from reciprocal space."""

    def __init__(self, experiment, sigma_b, sigma_m, n_sigma):
        """Initialise with models and parameters."""
        self.experiment = experiment
        self.sigma_b = sigma_b
        self.sigma_m = sigma_m
        self.n_sigma = n_sigma

    def with_given_intensity(self, N, In, Ba, Bb, Bc, Bd):
        """Generate reflections with a given intensity and background."""
        return self.with_individual_given_intensity(
            N,
            flex.int(N, In),
            flex.int(N, Ba),
            flex.int(N, Bb),
            flex.int(N, Bc),
            flex.int(N, Bd),
        )

    def with_random_intensity(self, N, Imax, Bamax, Bbmax, Bcmax, Bdmax):
        """Generate reflections with a random intensity and background."""
        if Imax == 0:
            In = flex.size_t(N).as_int()
        else:
            In = flex.random_size_t(N, Imax).as_int()
        if Bamax == 0:
            Ba = flex.size_t(N).as_int()
        else:
            Ba = flex.random_size_t(N, Bamax).as_int()
        if Bbmax == 0:
            Bb = flex.size_t(N).as_int()
        else:
            Bb = flex.random_size_t(N, Bbmax).as_int()
        if Bcmax == 0:
            Bc = flex.size_t(N).as_int()
        else:
            Bc = flex.random_size_t(N, Bcmax).as_int()
        if Bdmax == 0:
            Bd = flex.size_t(N).as_int()
        else:
            Bd = flex.random_size_t(N, Bdmax).as_int()
        return self.with_individual_given_intensity(N, In, Ba, Bb, Bc, Bd)

    def with_individual_given_intensity(self, N, In, Ba, Bb, Bc, Bd):
        """Generate reflections with given intensity and background."""
        from dials.algorithms.simulation import simulate_reciprocal_space_gaussian
        from dials.algorithms.simulation.generate_test_reflections import (
            random_background_plane2,
        )
        from dials.util.command_line import ProgressBar

        # Check the lengths
        assert N == len(In)
        assert N == len(Ba)
        assert N == len(Bb)
        assert N == len(Bc)
        assert N == len(Bd)

        # Generate some predictions
        refl = self.generate_predictions(N)

        # Calculate the signal
        progress = ProgressBar(title=f"Calculating signal for {len(refl)} reflections")
        s1 = refl["s1"]
        phi = refl["xyzcal.mm"].parts()[2]
        bbox = refl["bbox"]
        shoebox = refl["shoebox"]
        m = int(len(refl) / 100)
        I_exp = flex.double(len(refl), 0)
        for i in range(len(refl)):
            if In[i] > 0:
                data = shoebox[i].data.as_double()
                I_exp[i] = simulate_reciprocal_space_gaussian(
                    self.experiment.beam,
                    self.experiment.detector,
                    self.experiment.goniometer,
                    self.experiment.scan,
                    self.sigma_b,
                    self.sigma_m,
                    s1[i],
                    phi[i],
                    bbox[i],
                    In[i],
                    data,
                    shoebox[i].mask,
                )
                shoebox[i].data = data.as_float()
            if i % m == 0:
                progress.update(100.0 * float(i) / len(refl))
        progress.finished(f"Calculated signal impacts for {len(refl)} reflections")

        # Calculate the background
        progress = ProgressBar(
            title=f"Calculating background for {len(refl)} reflections"
        )
        for l in range(len(refl)):
            background = flex.float(flex.grid(shoebox[l].size()), 0.0)
            random_background_plane2(background, Ba[l], Bb[l], Bc[l], Bd[l])
            shoebox[l].data += background
            shoebox[l].background = background
            if l % m == 0:
                progress.update(100.0 * float(l) / len(refl))
            progress.update(100.0 * float(l) / len(refl))
        progress.finished(f"Calculated background for {len(refl)} reflections")

        ## Calculate the expected intensity by monte-carlo integration
        # progress = ProgressBar(title='Integrating expected signal for %d reflections' % len(refl))
        # s1 = refl['s1']
        # phi = refl['xyzcal.mm'].parts()[2]
        # bbox = refl['bbox']
        # shoebox = refl['shoebox']
        # I_exp = flex.double(len(refl), 0)
        # m = int(len(refl) / 100)
        # for i in range(len(refl)):
        # if In[i] > 0:
        # I_exp[i] = integrate_reciprocal_space_gaussian(
        # self.experiment.beam,
        # self.experiment.detector,
        # self.experiment.goniometer,
        # self.experiment.scan,
        # self.sigma_b,
        # self.sigma_m,
        # s1[i],
        # phi[i],
        # bbox[i],
        # 10000,
        # shoebox[i].mask) / 10000.0
        # if i % m == 0:
        # progress.update(100.0 * float(i) / len(refl))
        # progress.finished('Integrated expected signal impacts for %d reflections' % len(refl))

        # Save the expected intensity and background
        refl["intensity.sim"] = In
        refl["background.sim.a"] = Ba
        refl["background.sim.b"] = Bb
        refl["background.sim.c"] = Bc
        refl["background.sim.d"] = Bd
        refl["intensity.exp"] = I_exp

        # Return the reflections
        return refl

    def generate_predictions(self, N):
        """Generate some reflections."""
        from dials.algorithms import filtering
        from dials.algorithms.profile_model.gaussian_rs import MaskCalculator3D
        from dials.algorithms.profile_model.gaussian_rs import Model as ProfileModel
        from dials.algorithms.shoebox import MaskCode
        from dials.util.command_line import Command

        # Set the profile model
        self.experiment.profile = ProfileModel(
            None, self.n_sigma, self.sigma_b, self.sigma_m
        )

        # Generate a list of reflections
        refl = flex.reflection_table.from_predictions(self.experiment)
        refl["id"] = flex.int(len(refl), 0)

        # Filter by zeta
        zeta = 0.05
        Command.start(f"Filtering by zeta >= {zeta:f}")
        mask = filtering.by_zeta(
            self.experiment.goniometer, self.experiment.beam, refl["s1"], zeta
        )
        refl.del_selected(~mask)
        Command.end(f"Filtered {len(refl)} reflections by zeta >= {zeta:f}")

        # Compute the bounding box
        refl.compute_bbox([self.experiment])
        index = []
        image_size = self.experiment.detector[0].get_image_size()
        array_range = self.experiment.scan.get_array_range()
        bbox = refl["bbox"]
        for i in range(len(refl)):
            x0, x1, y0, y1, z0, z1 = bbox[i]
            if (
                x0 < 0
                or x1 > image_size[0]
                or y0 < 0
                or y1 > image_size[1]
                or z0 < array_range[0]
                or z1 > array_range[1]
            ):
                index.append(i)
        refl.del_selected(flex.size_t(index))

        # Sample if specified
        index = random.sample(range(len(refl)), N)
        refl = refl.select(flex.size_t(index))

        # Compute the bounding box
        # Create a load of shoeboxes
        Command.start(f"Creating shoeboxes for {len(refl)} reflections")
        refl["shoebox"] = flex.shoebox(refl["panel"], refl["bbox"])
        refl["shoebox"].allocate_with_value(MaskCode.Valid)
        Command.end(f"Created shoeboxes for {len(refl)} reflections")

        # Get the function object to mask the foreground
        Command.start(f"Masking Foreground for {len(refl)} reflections")
        mask_foreground = MaskCalculator3D(
            self.experiment.beam,
            self.experiment.detector,
            self.experiment.goniometer,
            self.experiment.scan,
            self.n_sigma * self.sigma_b,
            self.n_sigma * self.sigma_m,
        )

        # Mask the foreground
        mask_foreground(
            refl["shoebox"], refl["s1"], refl["xyzcal.px"].parts()[2], refl["panel"]
        )
        Command.end(f"Masked foreground for {len(refl)} reflections")

        # Return the reflections
        return refl
