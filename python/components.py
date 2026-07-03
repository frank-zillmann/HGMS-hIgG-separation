"""Component system for the HGMS hIgG separation model."""

import fs3

# Surface-site density and specific surface area of the MNPs (used by reactions and solutions).
mnp_ns = 1.08e-5
mnp_specific_surface = 1e5


def build_component_system() -> fs3.ComponentSystem:
    return fs3.ComponentSystem(
        [
            fs3.Component("H₂O", molar_mass_kg_per_mol=18.01528e-3),
            fs3.Component("H⁺", charge=1, molar_mass_kg_per_mol=1.008e-3, truesdell_jones_alpha=4.78e-10, truesdell_jones_beta=0.24e-3),
            fs3.Component("OH⁻", charge=-1, molar_mass_kg_per_mol=17.007e-3, truesdell_jones_alpha=10.65e-10, truesdell_jones_beta=0.21e-3),
            fs3.Component("Na⁺", charge=1, molar_mass_kg_per_mol=22.990e-3, truesdell_jones_alpha=4.32e-10, truesdell_jones_beta=0.06e-3),
            fs3.Component("Cl⁻", charge=-1, molar_mass_kg_per_mol=35.453e-3, truesdell_jones_alpha=3.71e-10, truesdell_jones_beta=0.01e-3),
            fs3.Component("TrisH⁺", charge=1, molar_mass_kg_per_mol=121.14e-3, truesdell_jones_alpha=4.0e-10, truesdell_jones_beta=0.0e-3),
            fs3.Component("Tris", molar_mass_kg_per_mol=121.14e-3),
            fs3.Component("AcH", molar_mass_kg_per_mol=60.052e-3),
            fs3.Component("Ac⁻", charge=-1, molar_mass_kg_per_mol=59.044e-3, truesdell_jones_alpha=4.5e-10, truesdell_jones_beta=0.0e-3),
            fs3.Component("GlyH₂⁺", charge=1, molar_mass_kg_per_mol=75.067e-3),
            fs3.Component("GlyH", molar_mass_kg_per_mol=75.067e-3),
            fs3.Component("Gly⁻", charge=-1, molar_mass_kg_per_mol=74.059e-3, truesdell_jones_alpha=4.0e-10, truesdell_jones_beta=0.0e-3),
            fs3.Component("hIgG", molar_mass_kg_per_mol=150.0),
            fs3.Component("MNP-OH₂⁺", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=18.015e-3),
            fs3.Component("MNP-OH", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=17.007e-3),
            fs3.Component("MNP-O⁻", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=15.999e-3),
            fs3.Component("MNP-hIgG", type=fs3.ComponentType.MagneticNanoParticleGroup, molar_mass_kg_per_mol=150000e-3),
            fs3.Component("MNP1", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1039e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP2", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1209e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP3", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1406e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP4", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1635e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP5", type=fs3.ComponentType.MagneticNanoParticle, radius_m=1901e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP6", type=fs3.ComponentType.MagneticNanoParticle, radius_m=2211e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP7", type=fs3.ComponentType.MagneticNanoParticle, radius_m=2571e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP8", type=fs3.ComponentType.MagneticNanoParticle, radius_m=2990e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP9", type=fs3.ComponentType.MagneticNanoParticle, radius_m=3477e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
            fs3.Component("MNP10", type=fs3.ComponentType.MagneticNanoParticle, radius_m=4043e-9, density_kg_per_m3=5170, magnetic_saturation_A_per_m=3.5e5),
        ]
    )
