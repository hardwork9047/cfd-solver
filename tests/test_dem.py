"""
Unit tests for the dem (Discrete Element Method) package.
"""

import numpy as np
import pytest

from dem import ParticleSystem


class TestParticleSystemInit:
    """Tests for ParticleSystem initialization."""

    def test_initialization(self):
        """Test that fields are initialized with correct shapes."""
        ps = ParticleSystem(n_particles=10, domain_size=(5.0, 5.0), particle_radius=0.1)

        assert ps.n_particles == 10
        assert ps.width == 5.0
        assert ps.height == 5.0
        assert ps.radius == 0.1
        assert ps.positions.shape == (10, 2)
        assert ps.velocities.shape == (10, 2)
        assert ps.time == 0.0
        assert ps.iteration == 0

    def test_initial_velocities_zero(self):
        """Particles start with zero velocity."""
        ps = ParticleSystem(n_particles=5)

        np.testing.assert_array_equal(ps.velocities, 0.0)

    def test_particle_mass(self):
        """Mass = density * pi * r² (2D circular cross-section)."""
        r = 0.2
        rho = 2500.0
        ps = ParticleSystem(n_particles=1, particle_radius=r, particle_density=rho)

        expected_mass = rho * np.pi * r**2
        assert np.isclose(ps.mass, expected_mass)

    def test_radii_and_masses_arrays(self):
        """radii and masses arrays must match the scalar radius and mass."""
        ps = ParticleSystem(n_particles=5, particle_radius=0.3, particle_density=1000.0)

        np.testing.assert_allclose(ps.radii, 0.3)
        np.testing.assert_allclose(ps.masses, ps.mass)

    def test_positions_within_domain(self):
        """Initial positions must lie inside the domain (accounting for radius)."""
        ps = ParticleSystem(n_particles=20, domain_size=(10.0, 10.0), particle_radius=0.1)

        assert np.all(ps.positions[:, 0] >= 0.0)
        assert np.all(ps.positions[:, 0] <= 10.0)
        assert np.all(ps.positions[:, 1] >= 0.0)
        assert np.all(ps.positions[:, 1] <= 10.0)


class TestGravity:
    """Tests for gravitational force."""

    def test_gravity_direction(self):
        """Gravity must act in the -y direction only."""
        ps = ParticleSystem(n_particles=3, gravity=9.81)

        # Place particles away from walls so only gravity acts
        ps.positions[:] = np.array([[5.0, 5.0], [5.0, 5.0], [5.0, 5.0]])
        ps.velocities[:] = 0.0

        forces = ps.compute_forces()

        # x-component of gravity must be zero
        np.testing.assert_allclose(forces[:, 0], 0.0, atol=1e-10)
        # y-component must be negative (downward)
        assert np.all(forces[:, 1] < 0.0)

    def test_gravity_magnitude(self):
        """Gravity force magnitude = mass * g."""
        g = 9.81
        ps = ParticleSystem(n_particles=2, gravity=g)

        ps.positions[:] = np.array([[5.0, 5.0], [5.0, 5.0]])
        ps.velocities[:] = 0.0

        forces = ps.compute_forces()

        expected_fy = -ps.mass * g
        np.testing.assert_allclose(forces[:, 1], expected_fy, rtol=1e-10)

    def test_zero_gravity(self):
        """With gravity=0 and no contacts, forces must be zero."""
        ps = ParticleSystem(n_particles=2, gravity=0.0)

        ps.positions[:] = np.array([[2.0, 5.0], [8.0, 5.0]])
        ps.velocities[:] = 0.0

        forces = ps.compute_forces()

        np.testing.assert_allclose(forces, 0.0, atol=1e-10)


class TestStep:
    """Tests for a single time step."""

    def test_time_advances(self):
        """Simulation time and iteration counter must advance after a step."""
        ps = ParticleSystem(n_particles=5, dt=1e-4)

        ps.step()

        assert np.isclose(ps.time, 1e-4)
        assert ps.iteration == 1

    def test_positions_change_under_gravity(self):
        """Particles fall under gravity — y-positions decrease over time."""
        ps = ParticleSystem(n_particles=1, gravity=9.81, dt=1e-3)

        ps.positions[:] = [[5.0, 5.0]]
        ps.velocities[:] = 0.0

        y_initial = ps.positions[0, 1]
        for _ in range(50):
            ps.step()

        assert ps.positions[0, 1] < y_initial, "Particle should fall under gravity."

    def test_no_nan_or_inf_after_step(self):
        """Fields must not contain NaN or Inf after a step."""
        ps = ParticleSystem(n_particles=10, dt=1e-4)

        ps.step()

        assert not np.any(np.isnan(ps.positions))
        assert not np.any(np.isnan(ps.velocities))
        assert not np.any(np.isinf(ps.positions))
        assert not np.any(np.isinf(ps.velocities))

    def test_particles_stay_in_domain(self):
        """After many steps, all particles must remain inside the domain."""
        ps = ParticleSystem(
            n_particles=5,
            domain_size=(5.0, 5.0),
            particle_radius=0.1,
            gravity=9.81,
            dt=1e-4,
        )

        for _ in range(200):
            ps.step()

        assert np.all(ps.positions[:, 0] >= 0.0)
        assert np.all(ps.positions[:, 0] <= ps.width)
        assert np.all(ps.positions[:, 1] >= 0.0)
        assert np.all(ps.positions[:, 1] <= ps.height)


class TestContactForce:
    """Tests for particle-particle contact forces."""

    def test_contact_force_is_repulsive(self):
        """Overlapping particles must push each other apart."""
        ps = ParticleSystem(n_particles=2, particle_radius=0.5, gravity=0.0)

        # Place particles so they overlap
        ps.positions[:] = np.array([[5.0, 5.0], [5.5, 5.0]])
        ps.velocities[:] = 0.0

        forces = ps.compute_forces()

        # Particle 0 should be pushed left (negative x), particle 1 pushed right
        assert forces[0, 0] < 0.0
        assert forces[1, 0] > 0.0

    def test_no_contact_no_force(self):
        """Well-separated particles must not exert contact forces on each other."""
        ps = ParticleSystem(n_particles=2, particle_radius=0.1, gravity=0.0)

        ps.positions[:] = np.array([[1.0, 5.0], [9.0, 5.0]])
        ps.velocities[:] = 0.0

        forces = ps.compute_forces()

        # With no contact and no gravity, all forces must be zero
        np.testing.assert_allclose(forces, 0.0, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
